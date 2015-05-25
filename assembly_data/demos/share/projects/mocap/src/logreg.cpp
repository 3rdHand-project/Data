#include "logreg.h"
#include <Optim/optimization.h>

LogReg::LogReg() { setDefaultParams(); }
LogReg::~LogReg() {};

void LogReg::setDefaultParams() {
  String thresh = MT::getParameter<String>("thresh");
  String optim = MT::getParameter<String>("optim");

  params.set<double>("cclip", 0.);
  params.set<double>("thresh", .5);
  params.set<String>("optim", optim);
}

void LogReg::update_p(DataLogReg &data, const arr &_theta) {
  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");
  if(data.last_theta_p.N != 0 && data.last_theta_p == theta)
    return;
  data.last_theta_p = theta;

  uint ndata = data.ndata;
  const arr &phi = data.phi;
  arr &phiTtheta = data.phiTtheta;
  arr &ephiTtheta = data.ephiTtheta;
  arr &logZ = data.logZ;
  arr &p = data.p;
  arr &p2d = data.p2d;

  ephiTtheta.resize(ndata);
  logZ.resize(ndata);
  p.resize(ndata);
  p2d.resize(ndata, 2);

  phiTtheta = phi * theta;

  double cclip = *params.get<double>("cclip");
  if(cclip > 0)
    for(uint i = 0; i < ndata; i++)
      ephiTtheta(i) = exp(clip(phiTtheta(i), -cclip, cclip));
  else
    for(uint i = 0; i < ndata; i++)
      ephiTtheta(i) = exp(phiTtheta(i));

  for(uint i = 0; i < ndata; i++) {
    logZ(i) = log1p(ephiTtheta(i));
    p(i) = ephiTtheta(i) / (1+ephiTtheta(i));
    p2d[i]() = ARR(1-p(i), p(i));
  }
}

Model::Data *LogReg::newData() {
  return new DataLogReg();
}

void LogReg::setFeatures(Data &data, const arr &d) {
  data.phi = d;
  params.set<uint>("nfeats", data.phi.d1);
}

double LogReg::negloglikelihood(const arr &x, Data &_data, const arr &_theta) {
  return 0./0.;
}

double LogReg::negloglikelihood(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  uint ndata = data.ndata;
  const arr &target = data.target;
  const arr &phiTtheta = data.phiTtheta;
  const arr &logZ = data.logZ;

  double nll = 0;
  for(uint i = 0; i < ndata; i++)
    nll += logZ(i) + (target(i)? -phiTtheta(i): 0);
  return nll;
}

arr LogReg::jacobian(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  const arr &target = data.target;
  const arr &phi = data.phi;
  const arr &p = data.p;

  return ~phi * (p - target);
}

arr LogReg::hessian(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  uint nfeats = data.nfeats;
  const arr &phi = data.phi;
  const arr &p = data.p;
  arr &p1p = data.p1p;

  p1p = p % (1. - p);

  arr H(nfeats, nfeats);
  H.setZero();
  for(uint i = 0; i < nfeats; i++)
    for(uint j = 0; j < nfeats; j++)
      H(i, j) = scalarProduct(phi[i] % p1p, phi[j]);
  return H;
}

arr LogReg::sample(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  double nframes = data.ndata;
  const arr &p = data.p;

  arr sample(nframes);
  for(uint t = 0; t < nframes; t++)
    sample(t) = (rnd.uni() <= p(t));
  return sample;
}

arr LogReg::maxlikelihood(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  double ndata = data.ndata;
  const arr &p = data.p;

  double thresh = *params.get<double>("thresh");

  arr ml(ndata);
  for(uint i = 1; i < ndata; i++)
    ml(i) = (thresh <= p(i));
  return ml;
}

arr LogReg::marginal(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);
  return data.p2d;
}

arr LogReg::marginal2(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  uint nframes = data.ndata;
  const arr &p2d = data.p2d;

  arr pp(nframes, 2, 2);
  for(uint i = 0; i < 2; i++)
    pp[0][i]() = p2d(0, i);
  for(uint t = 1; t < nframes; t++)
    pp[t]() = p2d[t] ^ p2d[t-1];
  return pp;
}

double LogReg::nseq_expected(Data &_data, const arr &_theta) {
  DataLogReg &data = static_cast<DataLogReg&>(_data);
  update_p(data, _theta);

  double nframes = data.ndata;
  const arr &p = data.p;

  double nseq = p(0);
  for(uint t = 1; t < nframes; t++)
    nseq += p(t) * (1-p(t-1));

  return nseq;
}

