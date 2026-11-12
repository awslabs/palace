// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "failurereport.hpp"

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
  bool active = false, root = false;
  FailureReport::Writer writer;
  double grace = 0.0;
  std::terminate_handler prev_terminate = nullptr;
  std::unique_ptr<TailBuffer> err_buf;
  std::streambuf *prev_err = nullptr;
#if defined(MFEM_USE_EXCEPTIONS)
  mfem::ErrorAction prev_action = mfem::MFEM_ERROR_THROW;
#endif
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
  if (!state.root)
  {
    // No collective is safe, as the other processes may not fail: give the root the time
    // to write its report, should it fail too (it then aborts all processes).
    std::this_thread::sleep_for(std::chrono::duration<double>(state.grace));
    return false;
  }
  std::string msg(message);
  if (msg.empty())
  {
    msg = state.err_buf->LastMessage(message_age);
  }
  msg.erase(msg.find_last_not_of(" \t\r\n") + 1);
  const nlohmann::json error = {{"Message", msg.empty() ? std::string(fallback) : msg},
                                {"Phase", Timer::names[BlockTimer::Current()]}};
  try
  {
    state.writer(error);
  }
  catch (...)
  {
    // The report is best effort.
    return false;
  }
  return true;
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

FailureReport::FailureReport(MPI_Comm comm, Writer writer, double grace)
{
  MFEM_VERIFY(!state.active, "Only one FailureReport may be active at a time!");
  state.root = Mpi::Root(comm);
  state.writer = std::move(writer);
  state.grace = grace;
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
}

bool ReportFailure(std::string_view message, std::string_view fallback)
{
  return Report(message, fallback);
}

}  // namespace palace::perf
