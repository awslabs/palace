#define _GNU_SOURCE
#include <dlfcn.h>
#include <execinfo.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef int cudaError_t;
typedef int CUresult;
typedef uintptr_t CUdeviceptr;
typedef void *cudaStream_t;

typedef struct Rec {
  void *ptr;
  size_t size;
  const char *api;
  int nbt;
  void *bt[12];
  struct Rec *next;
} Rec;

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static Rec *head = NULL;
static unsigned long long live_bytes = 0, peak_live_bytes = 0;
static unsigned long long alloc_calls = 0, free_calls = 0, failed_alloc_calls = 0;
static __thread int in_hook = 0;

static void *nextsym(const char *name) {
  void *p = dlsym(RTLD_NEXT, name);
  if (!p) {
    fprintf(stderr, "[cuda_memtrace] missing %s: %s\n", name, dlerror());
    abort();
  }
  return p;
}

static void add_rec(void *ptr, size_t size, const char *api) {
  if (!ptr) return;
  Rec *r = (Rec *)malloc(sizeof(Rec));
  if (!r) return;
  r->ptr = ptr;
  r->size = size;
  r->api = api;
  r->nbt = backtrace(r->bt, 12);
  pthread_mutex_lock(&lock);
  r->next = head;
  head = r;
  live_bytes += size;
  if (live_bytes > peak_live_bytes) peak_live_bytes = live_bytes;
  alloc_calls++;
  pthread_mutex_unlock(&lock);
}

static void remove_rec(void *ptr) {
  if (!ptr) return;
  pthread_mutex_lock(&lock);
  Rec **pp = &head;
  while (*pp) {
    Rec *r = *pp;
    if (r->ptr == ptr) {
      *pp = r->next;
      live_bytes -= r->size;
      free_calls++;
      pthread_mutex_unlock(&lock);
      free(r);
      return;
    }
    pp = &r->next;
  }
  free_calls++;
  pthread_mutex_unlock(&lock);
}

static void print_meminfo(void) {
  typedef cudaError_t (*F)(size_t *, size_t *);
  static F fn = NULL;
  if (!fn) fn = (F)dlsym(RTLD_NEXT, "cudaMemGetInfo");
  if (fn) {
    size_t free_b = 0, total_b = 0;
    if (fn(&free_b, &total_b) == 0) {
      fprintf(stderr, "[cuda_memtrace] cudaMemGetInfo free=%zu total=%zu used=%zu\n",
              free_b, total_b, total_b - free_b);
    }
  }
}

static void dump_live_locked(size_t limit) {
  fprintf(stderr,
          "[cuda_memtrace] summary live_bytes=%llu peak_live_bytes=%llu alloc_calls=%llu free_calls=%llu failed_alloc_calls=%llu\n",
          live_bytes, peak_live_bytes, alloc_calls, free_calls, failed_alloc_calls);
  size_t n = 0;
  for (Rec *r = head; r && n < limit; r = r->next, n++) {
    fprintf(stderr, "[cuda_memtrace] live[%zu] ptr=%p size=%zu api=%s\n", n, r->ptr,
            r->size, r->api);
    backtrace_symbols_fd(r->bt, r->nbt, fileno(stderr));
  }
}

static void report_failure(const char *api, size_t size, int rc) {
  pthread_mutex_lock(&lock);
  failed_alloc_calls++;
  fprintf(stderr, "[cuda_memtrace] FAILED %s size=%zu rc=%d live_bytes=%llu peak_live_bytes=%llu\n",
          api, size, rc, live_bytes, peak_live_bytes);
  dump_live_locked(40);
  pthread_mutex_unlock(&lock);
  print_meminfo();
  fflush(stderr);
}

cudaError_t cudaMalloc(void **ptr, size_t size) {
  typedef cudaError_t (*F)(void **, size_t);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cudaMalloc");
  cudaError_t rc = fn(ptr, size);
  if (in_hook) return rc;
  in_hook = 1;
  if (rc == 0) add_rec(ptr ? *ptr : NULL, size, "cudaMalloc");
  else report_failure("cudaMalloc", size, rc);
  in_hook = 0;
  return rc;
}

cudaError_t cudaFree(void *ptr) {
  typedef cudaError_t (*F)(void *);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cudaFree");
  if (!in_hook) {
    in_hook = 1;
    remove_rec(ptr);
    in_hook = 0;
  }
  return fn(ptr);
}

cudaError_t cudaMallocAsync(void **ptr, size_t size, cudaStream_t stream) {
  typedef cudaError_t (*F)(void **, size_t, cudaStream_t);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cudaMallocAsync");
  cudaError_t rc = fn(ptr, size, stream);
  if (in_hook) return rc;
  in_hook = 1;
  if (rc == 0) add_rec(ptr ? *ptr : NULL, size, "cudaMallocAsync");
  else report_failure("cudaMallocAsync", size, rc);
  in_hook = 0;
  return rc;
}

cudaError_t cudaFreeAsync(void *ptr, cudaStream_t stream) {
  typedef cudaError_t (*F)(void *, cudaStream_t);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cudaFreeAsync");
  if (!in_hook) {
    in_hook = 1;
    remove_rec(ptr);
    in_hook = 0;
  }
  return fn(ptr, stream);
}

CUresult cuMemAlloc_v2(CUdeviceptr *dptr, size_t bytesize) {
  typedef CUresult (*F)(CUdeviceptr *, size_t);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cuMemAlloc_v2");
  CUresult rc = fn(dptr, bytesize);
  if (in_hook) return rc;
  in_hook = 1;
  if (rc == 0) add_rec((void *)(uintptr_t)(dptr ? *dptr : 0), bytesize, "cuMemAlloc_v2");
  else report_failure("cuMemAlloc_v2", bytesize, rc);
  in_hook = 0;
  return rc;
}

CUresult cuMemAlloc(CUdeviceptr *dptr, size_t bytesize) {
  typedef CUresult (*F)(CUdeviceptr *, size_t);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cuMemAlloc");
  CUresult rc = fn(dptr, bytesize);
  if (in_hook) return rc;
  in_hook = 1;
  if (rc == 0) add_rec((void *)(uintptr_t)(dptr ? *dptr : 0), bytesize, "cuMemAlloc");
  else report_failure("cuMemAlloc", bytesize, rc);
  in_hook = 0;
  return rc;
}

CUresult cuMemFree_v2(CUdeviceptr dptr) {
  typedef CUresult (*F)(CUdeviceptr);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cuMemFree_v2");
  if (!in_hook) {
    in_hook = 1;
    remove_rec((void *)(uintptr_t)dptr);
    in_hook = 0;
  }
  return fn(dptr);
}

CUresult cuMemFree(CUdeviceptr dptr) {
  typedef CUresult (*F)(CUdeviceptr);
  static F fn = NULL;
  if (!fn) fn = (F)nextsym("cuMemFree");
  if (!in_hook) {
    in_hook = 1;
    remove_rec((void *)(uintptr_t)dptr);
    in_hook = 0;
  }
  return fn(dptr);
}

__attribute__((destructor)) static void report(void) {
  if (in_hook) return;
  in_hook = 1;
  pthread_mutex_lock(&lock);
  dump_live_locked(20);
  pthread_mutex_unlock(&lock);
  print_meminfo();
  fflush(stderr);
  in_hook = 0;
}
