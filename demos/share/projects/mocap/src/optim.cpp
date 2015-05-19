#include <Core/array.h>
#include <Core/util.h>
#include "optim.h"

double LineSearch::search(const arr &x, const arr &d, double alpha) {
  double fx = objf.f(x);
  arr J = objf.jacobian(x);
  double J_T_d = scalarProduct(J, d);
  for(;;) {
    // cout << "alpha: " << alpha << endl;
    if(objf.f(x + alpha * d) <= fx + rho_ls * alpha * J_T_d)
      break;
    alpha *= rho_alpha_minus;
  }
  return alpha;
}

void GradientDescent::loopUntil(arr &x, double end_condition) {
  arr J, d, delta;
  uint count = 0;

  double fx, alpha = 1;

  LineSearch ls;
  ls.objf = objf;
  for(uint t = 0;; t++) {
    fx = objf.f(x);
    J = objf.jacobian(x);

    d = -J/length(J);
    alpha = ls.search(x, d, alpha);
    delta = alpha * d;

    cout << "==========================" << endl;
    cout << "LOOP " << t << endl;
    cout << "count " << count << endl;
    cout << "x: " << x << endl;
    cout << "f(x): " << fx << endl;
    cout << "J(x): " << J << endl;
    cout << "d: " << d << endl;
    cout << "alpha: " << alpha << endl;
    cout << "delta: " << delta << endl;

    x += delta;
    // alpha = 1;
    alpha *= rho_alpha_plus;

    if(absMax(delta) < end_condition)
      break;
    // if(absMax(delta) >= end_condition)
    //   count = 0;
    // else {
    //   count++;
    //   if(count == 10)
    //     break;
    // }
  }
}

void BFGS::loopUntil(arr &x, double end_condition) {
  arr I, newJ, d, delta, newx, y, tmpH, tmpHT;
  double step, delta_T_y;

  uint n = x.N;

  I = eye(n);
  invH = I;

  LineSearch ls;
  ls.objf = objf;

  if(objf.f) objf.f(x);
  arr J = objf.jacobian(x);
  for(uint t = 0;; t++) {
    // cout << "LOOP " << t << endl;
    // cout << "x: " << x << endl;
    // cout << "f(x): " << objf.f(x) << endl;
    // cout << "J: " << J << endl;
    d = - invH * J;
    d /= length(d);

    // cout << "invH: " << invH << endl;
    // cout << "d: " << d << endl;

    step = ls.search(x, d);
    delta = step * d;
    // cout << "step: " << step << endl;
    // cout << "delta: " << delta << endl;
    newx = x + delta;
    if(objf.f) objf.f(newx);
    newJ = objf.jacobian(newx);

    y = newJ - J;
    x = newx;
    J = newJ;

    delta_T_y = scalarProduct(delta, y);
    tmpH = (I - (y ^ delta)) / delta_T_y;
    transpose(tmpHT, tmpH);
    invH = tmpHT * invH * tmpH + (delta ^ delta) / delta_T_y;

    // cout << "********" << endl;
    // cout << "x: " << x << endl;
    // cout << "delta: " << delta << endl;
    if(absMax(delta) < end_condition)
      break;
  }
}

void Newton::loopUntil(arr &x, double end_condition) {
  arr J, H, E;
  arr d, delta;
  double alpha;

  LineSearch ls;
  ls.objf = objf;

  E = diag(1, x.N);
  alpha = 1;
  for(uint t = 0;; t++) {
    if(objf.f) objf.f(x);
    J = objf.jacobian(x);
    H = objf.hessian(x);
    d = -inverse_SymPosDef(H + E) * J;
    alpha = ls.search(x, d, alpha);
    delta = alpha * d;
    alpha = MT::MIN(alpha * rho_alpha_plus, 1);

    cout << "LOOP " << t << endl;
    cout << "x: " << x << endl;
    if(objf.f) cout << "f(x): " << objf.f(x) << endl;
    cout << "J: " << J << endl;
    cout << "H: " << H << endl;
    cout << "d: " << d << endl;

    x += delta;
    if(absMax(delta) < end_condition)
      break;
  }
}

