#pragma once
#include <Core/array.h>
#include "model.h"

struct LogReg: Model {
  struct DataLogReg: Data {
    arr last_theta_p;
    arr last_theta_p1p;

    arr phiTtheta, ephiTtheta;
    arr logZ;
    arr p, p2d, p1p;
  };

  LogReg();
  ~LogReg();

  void setDefaultParams();

  void update_p(DataLogReg &data, const arr &theta = NoArr);

  Data *newData();
  void setFeatures(Data &data, const arr &d);
  double negloglikelihood(const arr &x, Data &data, const arr &theta = NoArr);
  double negloglikelihood(Data &data, const arr &theta = NoArr);
  arr jacobian(Data &data, const arr &theta = NoArr);
  arr hessian(Data &data, const arr &theta = NoArr);
 
  arr sample(Data &data, const arr &theta = NoArr);
  arr maxlikelihood(Data &data, const arr &theta = NoArr);
  arr marginal(Data &data, const arr &theta = NoArr);
  arr marginal2(Data &data, const arr &theta = NoArr);
  double nseq_expected(Data &data, const arr &theta = NoArr);
};
