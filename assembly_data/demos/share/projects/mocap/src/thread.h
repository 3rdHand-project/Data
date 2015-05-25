#pragma once
#include <Core/array.h>
#include <Core/keyValueGraph.h>
#include <Core/thread.h>

struct Job;
typedef MT::Array<Job*> JobL;

struct Worker;
typedef MT::Array<Worker*> WorkerL;

struct Job {
  std::function<void(void)> job;

  // uint state;
  // ConditionVariable cond;

  Job();
  ~Job();

  void operator()();
  // void wait();
};

struct Pool {
  ConditionVariable state;
  Mutex mutex;

  WorkerL workerl;
  JobL jobl;

  Pool();
  ~Pool();

  void set_pool_size(uint s);
  void init_pool();
  void destroy_pool();
  bool pool_on();
  void append_job(Job *job);
  void append_job(const JobL &jobl);
  Job *pop_job();
  void wait();
};

struct Worker: Thread {
  Pool &pool;

  Worker(Pool &_pool);
  ~Worker();

  Worker(Worker&&) = delete;
  Worker(Worker const&) = delete;
  Worker& operator=(Worker&&) = delete;
  Worker& operator=(Worker const&) = delete;

  void open();
  void close();
  void step();
};
