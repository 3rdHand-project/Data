#pragma once
#include <Core/array.h>
#include <Core/thread.h>
#include <Core/keyValueGraph.h>
#include <Optim/optimization.h>
#include "optim.h"
#include "thread.h"

struct Model: Parametric {
  String param_hash;
  // TODO what is this?
  ofstream file;
  Pool pool;

  struct Data {
    uint ndata, nfeats;
    arr target;
    arr phi;
  };
  typedef MT::Array<Data*> DataL;

  Model();
  virtual ~Model();

  virtual Data *newData() = 0;
  virtual void setFeatures(Data &data, const arr &d) = 0;
  virtual double negloglikelihood(const arr &x, Data &data, const arr &theta = NoArr) = 0;;
  virtual double negloglikelihood(Data &data, const arr &theta = NoArr) = 0;
  virtual arr jacobian(Data &data, const arr &theta = NoArr) = 0;
  virtual arr hessian(Data &data, const arr &theta = NoArr) = 0;

  virtual arr sample(Data &data, const arr &theta = NoArr) = 0;
  virtual arr maxlikelihood(Data &data, const arr &theta = NoArr) = 0;
  virtual arr marginal(Data &data, const arr &theta = NoArr) = 0;
  virtual arr marginal2(Data &data, const arr &theta = NoArr) = 0;
  virtual double nseq_expected(Data &data, const arr &theta = NoArr) = 0;

  double nseq_maxlikelihood(Data &data, const arr &theta = NoArr);
  double entropy(Data &data, const arr &theta = NoArr);

  double negloglikelihood(const DataL &datal, const arr &theta);
  arr jacobian(const DataL &data, const arr &theta);
  arr hessian(const DataL &data, const arr &theta);

  void saveParams();
  bool loadParams();
  void train(const DataL &datal);

  ObjectiveFunction makeObj(const DataL &datal);
};

