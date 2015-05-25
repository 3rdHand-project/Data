#pragma once

#include <Core/array.h>

struct ObjectiveFunction {
  uint arity;
  std::function<double(arr)> f;
  std::function<arr(arr)> jacobian, hessian;

  ObjectiveFunction() {};
  ObjectiveFunction(uint _arity): arity(_arity) {};
  ObjectiveFunction(const ObjectiveFunction &objf) {
    arity = objf.arity;
    f = objf.f;
    jacobian = objf.jacobian;
    hessian = objf.hessian;
  };
  ~ObjectiveFunction() {};

#define LIMIT 10.
  void testPoints(uint n = 10) {
    arr x, J;
    double fx;

    for(uint i = 0; i < n; i++) {
      x = LIMIT * (rand(TUP(arity)) - .5);
      fx = f(x);
      J = jacobian(x);
      
      cout << "=================" << endl;
      cout << "x: " << x << endl;
      cout << "fx: " << fx << endl;
      cout << "J: " << J << endl;
    }
  }

  bool convexCheck(uint n = 100) {
    arr x, y, dxy;
    double fx, fy, dfxy, ifxy, fixy;

    uint n_wrong = 0;
    for(uint i = 0; i < n; i++) {
      cout << "Checking Convexity. i: " << i << endl;
      x = LIMIT * (rand(TUP(arity)) - .5);
      y = LIMIT * (rand(TUP(arity)) - .5);
      fx = f(x);
      fy = f(y);
      dxy = y - x;
      dfxy = fy - fx;
      for(double t = .2; t < .8; t += .1) {
        ifxy = fx + t * dfxy;
        fixy = f(x + t * dxy);
        if(isnan(ifxy) || isnan(ifxy) || ifxy<fixy) {
          cout << "==============================" << endl;
          cout << "x: " << x << endl;
          cout << "y: " << y << endl;
          cout << "fx: " << fx << endl;
          cout << "fy: " << fy << endl;
          cout << "t: " << t << endl;
          cout << "ifxy: " << ifxy << endl;
          cout << "fixy: " << fixy << endl;
          n_wrong++;
        }
      }
    }

    cout << "Wrong convexity: " << n_wrong << endl;
    return n_wrong == 0;
  }

  bool gradientCheck(uint n = 100, double eps = 1e-4) {
    arr x, J, J_approx(arity);
    double fx;
    
    arr epsI = eps * eye(arity);
    double inv2Eps = 1. / (2. * eps);
    // cout << "inv2Eps: " << inv2Eps << endl;

    uint n_wrong = 0;
    for(uint i = 0; i < n; i++) {
      cout << "Checking Gradient. i: " << i << endl;
      x = LIMIT * (rand(TUP(arity)) - .5);
      for(uint d = 0; d < arity; d++)
        J_approx(d) = inv2Eps * (f(x + epsI[d]) - f(x - epsI[d]));
      fx = f(x);
      J = jacobian(x);
      if(J!=J || J_approx!=J_approx || absMax(J - J_approx) >= 1e-2) {
        cout << "===========================================" << endl;
        cout << "x: " << x << endl;
        cout << "f(x): " << fx << endl;
        cout << "J: " << J << endl;
        cout << "J_approx: " << J_approx << endl;
        cout << "J - J_approx: " << J - J_approx << endl;
        cout << "absMax J - J_approx: " << absMax(J - J_approx) << endl;
        cout << "absMax J - J_approx >= 1e-2: " << (absMax(J - J_approx) >= 1e-2) << endl;
        n_wrong++;
      }
    }

    cout << "Wrong gradients: " << n_wrong << endl;
    return n_wrong == 0;
  }

  bool hessianCheck(uint n = 100, double eps = 1e-4) {
    arr x, J, H, H_approx(arity, arity);
    double fx;
    
    arr epsI = eps * eye(arity);
    double inv4Eps2 = 1. / (4. * eps * eps);
    // cout << "inv4Eps2: " << inv4Eps2 << endl;

    uint n_wrong = 0;
    for(uint i = 0; i < n; i++) {
      cout << "Checking Hessian. i: " << i << endl;
      x = LIMIT * (rand(TUP(arity)) - .5);
      for(uint d1 = 0; d1 < arity; d1++)
        for(uint d2 = 0; d2 < arity; d2++)
          H_approx(d1, d2) = inv4Eps2 * 
            ( f(x + epsI[d1] + epsI[d2]) + f(x - epsI[d1] - epsI[d2])
            - f(x + epsI[d1] - epsI[d2]) - f(x - epsI[d1] + epsI[d2]));
      fx = f(x);
      J = jacobian(x);
      H = hessian(x);
      if(absMax(H - H_approx) >= 1e-2) {
        cout << "===========================================" << endl;
        cout << "x: " << x << endl;
        cout << "f(x): " << fx << endl;
        cout << "J: " << J << endl;
        // cout << "H: " << H << endl;
        // cout << "H_approx: " << H_approx << endl;
        cout << "H - H_approx: " << H - H_approx << endl;
        cout << "absMax H - H_approx: " << absMax(H - H_approx) << endl;
        cout << "absMax H - H_approx >= 1e-2: " << (absMax(H - H_approx) >= 1e-2) << endl;
        n_wrong++;
      }
    }

    cout << "Wrong hessians: " << n_wrong << endl;
    return n_wrong == 0;
  }

  void test(const arr &x) {
    if(f)
      cout << "f: " << f(x) << endl;
    if(jacobian)
      cout << "jacobian: " << jacobian(x) << endl;
    if(hessian)
      cout << "hessian: " << hessian(x) << endl;
  }
};

struct GradientDescent {
  ObjectiveFunction objf;
  double rho_alpha_plus;

  GradientDescent(): rho_alpha_plus(1.2) {}
  ~GradientDescent() {}

  void loopUntil(arr &x, double end_condition);
};

struct LineSearch {
  ObjectiveFunction objf;
  double alpha, rho_alpha_minus, rho_ls;

  LineSearch(): rho_alpha_minus(.5), rho_ls(.15) {}
  ~LineSearch() {}

  double search(const arr &x, const arr &d, double alpha = 10);
};

struct BFGS {
  ObjectiveFunction objf;
  arr invH;

  void loopUntil(arr &x, double end_condition);
};

struct Newton {
  double rho_alpha_plus;
  ObjectiveFunction objf;

  Newton(): rho_alpha_plus(1.2) {}
  ~Newton() {}

  void loopUntil(arr &x, double end_condition);
};

