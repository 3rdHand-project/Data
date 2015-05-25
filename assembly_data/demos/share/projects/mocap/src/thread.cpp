#include "thread.h"

// =============================================================================
// Job
//

// Job::Job(): state(0) {}
Job::Job() {}
Job::~Job() {}

void Job::operator()() {
  CHECK(job, "Job is empty!");
  job();
}

// =============================================================================
// Pool
// 

Pool::Pool() { workerl.resize(4); state.setValue(-1); }
Pool::~Pool() { destroy_pool(); }

void Pool::set_pool_size(uint s) {
  CHECK(!pool_on(), "Can't change pool size after init");
  workerl.resize(s);
}

void Pool::init_pool() {
  if(pool_on())
    return;

  state.setValue(0);

  CHECK(workerl.N, "Select a positive number of workers");
  for(uint i = 0; i < workerl.N; i++) {
    Worker *w = new Worker(*this);
    w->threadLoop();
    workerl(i) = w;
  }
  cout << "Pool created: " << workerl.N << " workers." << endl;
}

void Pool::destroy_pool() {
  if(pool_on()) {
    state.setValue(-1);
    for(Worker *w: workerl)
      w->threadClose();
  }
}

bool Pool::pool_on() {
  state.lock();
  int value = state.value;
  state.unlock();
  return value != -1;
}

void Pool::append_job(Job *job) {
  mutex.lock();
  jobl.append(job);
  mutex.unlock();

  state.incrementValue();
}

void Pool::append_job(const JobL &_jobl) {
  mutex.lock();
  jobl.append(_jobl);
  mutex.unlock();

  state.mutex.lock();
  state.value = state.value + _jobl.N;
  state.broadcast();
  state.mutex.unlock();
}

Job *Pool::pop_job() {
  mutex.lock();
  Job *job = jobl.N? jobl.popFirst(): nullptr;
  mutex.unlock();

  return job;
}

void Pool::wait() {
  state.waitForValueEq(0);
}

// =============================================================================
// Worker
//

Worker::Worker(Pool &_pool): Thread("Worker"), pool(_pool) {}
Worker::~Worker() {}

void Worker::open() { }
void Worker::close() { }
void Worker::step() {
  // pool.state.waitForValueNotEq(0);
  // pool.state.waitForSignal(1.);//ForValueNotEq(0);
  if(!pool.pool_on())
    return;
  Job *job = pool.pop_job();
  if(job == nullptr)
    return;

  (*job)();

  pool.state.lock();
  pool.state.value = pool.state.value-1;
  pool.state.broadcast();
  pool.state.unlock();
  delete job;
}


