// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "failurereport.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <exception>
#include <memory>
#include <mutex>
#include <streambuf>
#include <string>
#include <thread>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

namespace palace::perf
{

namespace
{

// Forwards to another stream buffer and keeps the last characters written, to recover the
// message MFEM prints just before it aborts.
class TailBuffer : public std::streambuf
{
  static constexpr std::size_t max_size = 4096;
  std::streambuf *dest;
  std::string tail;
  Timer::TimePoint last_write;
  bool written = false;
  mutable std::mutex mutex;

  void Keep(const char *s, std::size_t n)
  {
    std::lock_guard<std::mutex> lock(mutex);
    tail.append(s, n);
    if (tail.size() > 2 * max_size)
    {
      tail.erase(0, tail.size() - max_size);
    }
    last_write = Timer::Now();
    written = true;
  }

protected:
  int_type overflow(int_type c) override
  {
    if (traits_type::eq_int_type(c, traits_type::eof()))
    {
      return traits_type::not_eof(c);
    }
    const char ch = traits_type::to_char_type(c);
    Keep(&ch, 1);
    return dest ? dest->sputc(ch) : c;
  }

  std::streamsize xsputn(const char *s, std::streamsize n) override
  {
    Keep(s, static_cast<std::size_t>(n));
    return dest ? dest->sputn(s, n) : n;
  }

  int sync() override { return dest ? dest->pubsync() : 0; }

public:
  explicit TailBuffer(std::streambuf *dest) : dest(dest) {}

  // The last message, if written at most max_age seconds ago (otherwise it is not the one
  // of this failure): MFEM writes "\n\n<message>\n", optionally followed by a backtrace.
  std::string LastMessage(double max_age) const
  {
    std::lock_guard<std::mutex> lock(mutex);
    if (!written || Timer::Duration(Timer::Now() - last_write).count() > max_age)
    {
      return {};
    }
    std::string s = tail.substr(0, tail.rfind("Backtrace:"));
    s.erase(s.find_last_not_of(" \t\r\n") + 1);
    if (const auto pos = s.rfind("\n\n"); pos != std::string::npos)
    {
      s.erase(0, pos + 2);
    }
    s.erase(0, s.find_first_not_of(" \t\r\n"));
    return s;
  }
};

struct State
{
  bool active = false;
  MPI_Comm comm = MPI_COMM_NULL;  // Duplicate, for the agreement step and the report
  FailureReport::Writer writer;
  double timeout = 0.0;
  std::terminate_handler prev_terminate = nullptr;
  std::unique_ptr<TailBuffer> err_buf;
  std::streambuf *prev_err = nullptr;
#if defined(MFEM_USE_EXCEPTIONS)
  mfem::ErrorAction prev_action = mfem::MFEM_ERROR_THROW;
#endif

  // Each process posts exactly one agreement reduction while the FailureReport is alive:
  // {1, 1 on the root} when it fails, {0, 0} in the destructor otherwise. The sums are the
  // number of processes which failed and whether the root did.
  bool posted = false, completed = false;
  MPI_Request request = MPI_REQUEST_NULL;
  int failed[2] = {0, 0}, num_failed[2] = {0, 0};
};

State state;
std::atomic_flag reporting = ATOMIC_FLAG_INIT;
std::atomic<std::thread::id> reporter;

// An MFEM error message is the one of this failure if printed this recently (seconds).
constexpr double message_age = 2.0;

bool MpiActive()
{
  int initialized = 0, finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  return initialized && !finalized;
}

void PostAgreement(bool failed)
{
  state.failed[0] = failed;
  state.failed[1] = failed && Mpi::Root(state.comm);
  state.posted = true;
  MPI_Iallreduce(state.failed, state.num_failed, 2, MPI_INT, MPI_SUM, state.comm,
                 &state.request);
}

bool AllFailed()
{
  return state.completed && state.num_failed[0] == Mpi::Size(state.comm);
}

// Wait for the agreement step to complete, up to the deadline.
void WaitAgreement(Timer::TimePoint deadline)
{
  while (!state.completed)
  {
    int done = 0;
    MPI_Test(&state.request, &done, MPI_STATUS_IGNORE);
    state.completed = done;
    if (done || Timer::Now() > deadline)
    {
      break;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
}

bool Report(std::string_view message, std::string_view fallback)
{
  if (!state.active)
  {
    return false;
  }
  if (reporting.test_and_set())
  {
    // Another thread reports: wait for it to abort the process. The same thread (a failure
    // in the report itself) aborts directly.
    while (reporter.load() != std::this_thread::get_id())
    {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    return false;
  }
  reporter.store(std::this_thread::get_id());

  // MPI and the performance report belong to the main thread.
  int main_thread = 0;
  MPI_Is_thread_main(&main_thread);
  if (!main_thread)
  {
    return false;
  }
  std::string msg(message);
  if (msg.empty())
  {
    msg = state.err_buf->LastMessage(message_age);
  }
  msg.erase(msg.find_last_not_of(" \t\r\n") + 1);
  nlohmann::json error = {{"Message", msg.empty() ? std::string(fallback) : msg},
                          {"Phase", Timer::names[BlockTimer::Current()]}};

  // The failure may have happened within a verification probe: report as usual.
  ClearPause();

  // Agreement step: wait for all processes to fail, up to the timeout. A non-root process
  // which gives up waits for up to twice as long more before it aborts, unless the
  // agreement shows that the root did not fail (see FailureReport).
  const auto t0 = Timer::Now();
  if (!state.posted)
  {
    PostAgreement(true);
  }
  const auto timeout =
      std::chrono::duration_cast<Timer::Clock::duration>(Timer::Duration(state.timeout));
  WaitAgreement(t0 + timeout);
  if (!AllFailed() && !Mpi::Root(state.comm))
  {
    const auto deadline = t0 + 3 * timeout;
    WaitAgreement(deadline);
    if (state.completed && !AllFailed() && state.num_failed[1] > 0)
    {
      // The root writes its own report: give it the time to (it aborts all processes), at
      // least the timeout from now, as it may have failed just before the deadline.
      std::this_thread::sleep_until(std::max(deadline, Timer::Now() + timeout));
    }
  }
  const bool collective = AllFailed();
  const double wait = Timer::Duration(Timer::Now() - t0).count();
  detail::comm_wait += wait;
  AddTally(RegisterTally("Wait/FailureReport"), wait, 1, wait);
  error["Collective"] = collective;
  if (collective)
  {
    try
    {
      state.writer(state.comm, error);
    }
    catch (...)
    {
      // The report is best effort. The other processes may wait in a collective of the
      // report which this one left: abort all of them now.
      PMPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Barrier(state.comm);  // No process aborts before the root has written
  }
  else if (Mpi::Root(state.comm))
  {
    try
    {
      state.writer(MPI_COMM_SELF, error);
    }
    catch (...)
    {
    }
  }
  return collective;
}

[[noreturn]] void OnTerminate()
{
  std::string msg = "std::terminate called";
  if (auto e = std::current_exception())
  {
    try
    {
      std::rethrow_exception(e);
    }
    catch (const std::exception &ex)
    {
      msg = ex.what();
    }
    catch (...)
    {
      msg = "Unknown exception";
    }
  }
  // MPI_Abort needs MPI and, off the main thread, MPI_THREAD_MULTIPLE.
  int main_thread = 0, provided = MPI_THREAD_SINGLE;
  if (MpiActive())
  {
    MPI_Is_thread_main(&main_thread);
    MPI_Query_thread(&provided);
  }
  if (!MpiActive() || (!main_thread && provided < MPI_THREAD_MULTIPLE))
  {
    if (state.prev_terminate)
    {
      state.prev_terminate();
    }
    std::abort();
  }
  std::fprintf(stderr, "\nTerminating after an uncaught exception: %s\n", msg.c_str());
  std::fflush(stderr);
  Report(msg, {});
  PMPI_Abort(MPI_COMM_WORLD, 1);
  std::abort();
}

}  // namespace

FailureReport::FailureReport(MPI_Comm comm, Writer writer, double timeout)
{
  MFEM_VERIFY(!state.active, "Only one FailureReport may be active at a time!");
  MPI_Comm_dup(comm, &state.comm);
  state.writer = std::move(writer);
  state.timeout = timeout;
  state.posted = state.completed = false;
  state.request = MPI_REQUEST_NULL;
  state.err_buf = std::make_unique<TailBuffer>(mfem::err.rdbuf());
  state.prev_err = mfem::err.rdbuf(state.err_buf.get());
#if defined(MFEM_USE_EXCEPTIONS)
  state.prev_action = mfem::get_error_action();
  mfem::set_error_action(mfem::MFEM_ERROR_THROW);
#endif
  state.prev_terminate = std::set_terminate(&OnTerminate);
  reporter.store(std::thread::id());
  reporting.clear();
  state.active = true;
}

FailureReport::~FailureReport()
{
  state.active = false;
  std::set_terminate(state.prev_terminate);
#if defined(MFEM_USE_EXCEPTIONS)
  mfem::set_error_action(state.prev_action);
#endif
  mfem::err.rdbuf(state.prev_err);
  state.err_buf.reset();
  state.writer = nullptr;
  if (std::uncaught_exceptions() == 0)
  {
    // Also completes the agreement step of a process which reported a failure and did not
    // abort, so that the communicator can be freed.
    if (!state.posted)
    {
      PostAgreement(false);
    }
    MPI_Wait(&state.request, MPI_STATUS_IGNORE);
    MPI_Comm_free(&state.comm);
  }
  state.comm = MPI_COMM_NULL;
}

bool ReportFailure(std::string_view message)
{
  return Report(message, "Unknown error");
}

}  // namespace palace::perf

#if defined(PALACE_WITH_MPI_ABORT_HOOK)
// Interpose MPI_Abort (profiling interface) to write the failure report first. The symbol
// is weak so that a static tool library's own MPI_Abort takes precedence at link time (the
// dynamic linker ignores weakness, see FailureReport).
extern "C" __attribute__((weak)) int MPI_Abort(MPI_Comm comm, int errorcode)
{
  palace::perf::Report({}, fmt::format("MPI_Abort called with error code {}", errorcode));
  return PMPI_Abort(comm, errorcode);
}
#endif
