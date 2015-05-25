#pragma once
#include <Core/array.h>
#include "model.h"

struct LinChainCRF: Model {
  struct DataCRF: Data {
    arr last_theta_psi;
    arr last_theta_p;
    arr last_theta_pp;
    arr last_theta_ppp;

    arr phiTtheta;
    arr psi, psiT;
    arr a, b;
    arr p, pp, ppp;
    double sumloga, sumlogb, logz;
    arr dlogz;
  };

  LinChainCRF();
  ~LinChainCRF();

  void setDefaultParams();

  void update_psi(DataCRF &data, const arr &theta = NoArr);
  void update_p(DataCRF &data, const arr &theta = NoArr);
  void update_pp(DataCRF &data, const arr &theta = NoArr);
  void update_ppp(DataCRF &data, const arr &theta = NoArr);

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

