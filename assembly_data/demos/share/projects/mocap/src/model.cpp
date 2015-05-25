#include "model.h"
#include <Mocap/utilAB.h>

// #define UNTHREADED

Model::Model() {}
Model::~Model() {}

double Model::nseq_maxlikelihood(Data &data, const arr &_theta) {
  arr ml = maxlikelihood(data, _theta);

  double nframes = data.ndata;

  uint nseq = ml(0);
  for(uint t = 1; t < nframes; t++)
    if(ml(t) && !ml(t-1))
      nseq++;
  return nseq;
}

double Model::entropy(Data &data, const arr &theta) {
  // uint nsamples = 1000;
  // uint H = 0;
  // for(uint n = 0; n < nsamples; n++) {
  //   arr s = sample(data, theta);
  //   double nll = negloglikelihood(s, data, theta);
  //   H -= nll * log(nll);
  // }
  // return H;
  return 0./0.;
}

double Model::negloglikelihood(const DataL &datal, const arr &_theta) {
  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");

#ifdef UNTHREADED
  double nll = 0;
  for(Data *data: datal)
    nll += negloglikelihood(*data, theta);
  return nll;
#else
  uint njobs = datal.N;
  arr nll(njobs);
  for(uint i = 0; i < njobs; i++) {
    Job *job = new Job();
    job->job = [&, i]() { nll(i) = negloglikelihood(*datal(i), theta); };
    pool.append_job(job);
  }
  pool.wait();
  return sum(nll);
#endif
}

int i = 0;
arr Model::jacobian(const DataL &datal, const arr &_theta) {
  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");
  uint nfeats = theta.N;

#ifdef UNTHREADED
  arr J(nfeats);
  J.setZero();
  for(Data *data: datal)
    J += jacobian(*data, theta);
  return J;
#else
  uint njobs = datal.N;
  arr J(njobs, nfeats);
  for(uint i = 0; i < njobs; i++) {
    Job *job = new Job();
    job->job = [&, i]() { J[i]() = jacobian(*datal(i), theta); };
    pool.append_job(job);
  }
  pool.wait();
  return sum(J, 0);
#endif
}

arr Model::hessian(const DataL &datal, const arr &_theta) {
  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");
  uint nfeats = theta.N;

#ifdef UNTHREADED
  arr H(nfeats);
  H.setZero();
  for(Data *data: datal)
    H += hessian(*data, theta);
  return H;
#else
  uint njobs = datal.N;
  arr H(njobs, nfeats, nfeats);
  for(uint i = 0; i < njobs; i++) {
    Job *job = new Job();
    job->job = [&, i]() { H[i]() = hessian(*datal(i), theta); };
    pool.append_job(job);
  }
  pool.wait();
  return sum(H, 0).reshape(nfeats, nfeats);
#endif
}

void Model::saveParams() {
  // TODO go and accumulate this type of stuff
  String pfname;
  pfname << "z.params_" << hash::get("params");

  ofstream pfile;
  MT::open(pfile, pfname);
  arr theta = *params.get<arr>("theta");
  pfile << theta;
  pfile.close();
}

bool Model::loadParams() {
  String pfname;
  pfname << "z.params_" << hash::get("params");

  ifstream pfile;
  try {
    MT::open(pfile, pfname);
    arr theta;
    pfile >> theta;
    params.set<arr>("theta", theta);

    pfile.close();
    return true;
  }
  catch(const char *e) {
  }

  return false;
}

void Model::train(const DataL &datal) {
  watch::push();

  // String params = MT::getParameter<String>("params");
  // hash::append(params, );

#ifndef UNTHREADED
  pool.init_pool();
#endif

  uint nfeats = *params.get<uint>("nfeats");
  arr theta(nfeats);
  theta.setZero();

  ObjectiveFunction objf = makeObj(datal);

  String check_points = MT::getParameter<String>("check_points");
  if(check_points == "true") {
    MT::rnd.seed(0);
    objf.testPoints();
    cout << "Check_points: DONE" << endl;
  }

  String check_convex = MT::getParameter<String>("check_convex");
  if(check_convex == "true") {
    MT::rnd.seed(0);
    CHECK(objf.convexCheck(), "Convexity check failed");
    cout << "Check_convex: DONE" << endl;
  }

  String check_gradient = MT::getParameter<String>("check_gradient");
  if(check_gradient == "true") {
    MT::rnd.seed(0);
    CHECK(objf.gradientCheck(), "Gradient check failed");
    cout << "Check_gradient: DONE" << endl;
  }

  String check_hessian = MT::getParameter<String>("check_hessian");
  if(check_hessian == "true") {
    MT::rnd.seed(0);
    CHECK(objf.hessianCheck(), "Hessian check failed");
    cout << "Check_hessian: DONE" << endl;
  }

  double end_condition = MT::getParameter<double>("end_condition");
  String optim = MT::getParameter<String>("optim");
  if(optim == "gradient_descent") {
    GradientDescent gd;
    gd.objf = objf;
    gd.loopUntil(theta, end_condition);
    cout << "GD theta: " << theta << endl;
    cout << " - f: " << objf.f(theta) << endl;
  }
  else if(optim == "newton") {
    Newton newt;
    newt.objf = objf;
    newt.loopUntil(theta, end_condition);
    cout << "NEWT theta: " << theta << endl;
    cout << " - f: " << objf.f(theta) << endl;
  }
  else if(optim == "rprop") {
    // struct ObjF: ScalarFunction {
    //   ObjectiveFunction objf;
    //   ObjF(ObjectiveFunction _objf): objf(_objf) {};
    //   ~ObjF() {};
    //   double fs(arr &g, arr &H, const arr &x) {
    //     if(&g) {
    //       cout << "==================== JACOBIAN" << endl;
    //       g = objf.jacobian(x);
    //     }
    //     if(&H) {
    //       cout << "==================== HESSIAN" << endl;
    //       H = objf.hessian(x);
    //     }
    //     cout << "==================== F" << endl;
    //     return objf.f(x);
    //   }
    // };
    // ObjF f(objf);
    ScalarFunction f = [&objf] (arr &g, arr &H, const arr &x) {
      if(&g) g = objf.jacobian(x);
      if(&H) H = objf.hessian(x);
      return objf.f(x);
    };

    optRprop(theta, f, NOOPT);
    cout << "RPROP theta: " << theta << endl;
    cout << " - f: " << objf.f(theta) << endl;
  }
  else HALT("optimization method '" << optim << "' not implemented yet.")

  cout << "DONE:" << endl;
  cout << "  - theta: " << theta << endl;
  cout << "  - f(theta): " << objf.f(theta) << endl;
  cout << "  - J(theta): " << objf.jacobian(theta) << endl;

  params.set<arr>("theta", theta);
  saveParams();

#ifndef UNTHREADED
  pool.destroy_pool();
#endif

  cout << __PRETTY_FUNCTION__ << " done in " << watch::pop() << " sec." << endl;
}

ObjectiveFunction Model::makeObj(const DataL &datal) {
  ObjectiveFunction objf;

  objf.arity = datal(0)->nfeats;
  objf.f = [this, &datal](const arr &theta) {
    return this->negloglikelihood(datal, theta);
  };
  objf.jacobian = [this, &datal](const arr &theta) {
    return this->jacobian(datal, theta);
  };
  objf.hessian = [this, &datal](const arr &theta) {
    return this->hessian(datal, theta);
  };

  return objf;
}

