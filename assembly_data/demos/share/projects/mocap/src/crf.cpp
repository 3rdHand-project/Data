#include "crf.h"
#include <Optim/optimization.h>

LinChainCRF::LinChainCRF() { setDefaultParams(); }
LinChainCRF::~LinChainCRF() {};

void LinChainCRF::setDefaultParams() {
  String optim = MT::getParameter<String>("optim");

  params.set<double>("cclip", 0.);
  params.set<String>("optim", optim);
}

void LinChainCRF::update_psi(DataCRF &data, const arr &_theta) {
  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");

  if(data.last_theta_psi.N != 0 && data.last_theta_psi == theta)
    return;
  data.last_theta_psi = theta;

  uint nframes = data.ndata;
  const arr &phi = data.phi;
  arr &phiTtheta = data.phiTtheta;
  arr &psi = data.psi;
  arr &psiT = data.psiT;

  phiTtheta.resize(nframes, 2, 2);
  psi.resize(nframes, 2, 2);
  psiT.resize(nframes, 2, 2);

  double cclip = *params.get<double>("cclip");
  if(cclip > 0) {
    for(uint t = 0; t < nframes; t++) {
      for(uint i = 0; i < 2; i++) {
        for(uint j = 0; j < 2; j++) {
          phiTtheta(t, i, j) = scalarProduct(phi[t][i][j], theta);
          psiT(t, j, i) = psi(t, i, j) = exp(clip(phiTtheta(t, i, j), -cclip, cclip));
        }
      }
    }
  }
  else {
    for(uint t = 0; t < nframes; t++) {
      for(uint i = 0; i < 2; i++) {
        for(uint j = 0; j < 2; j++) {
          phiTtheta(t, i, j) = scalarProduct(phi[t][i][j], theta);
          psiT(t, j, i) = psi(t, i, j) = exp(phiTtheta(t, i, j));
        }
      }
    }
  }
  // for(uint t = 0; t < nframes; t++) {
  //   cout << "phi[" << t << "]: " << phi[t] << endl;
  //   cout << "phiTtheta[" << t << "]: " << phiTtheta[t] << endl;
  // }
}

void LinChainCRF::update_p(DataCRF &data, const arr &_theta) {
  update_psi(data, _theta);

  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");
  if(data.last_theta_p.N != 0 && data.last_theta_p == theta)
    return;
  data.last_theta_p = theta;

  uint nframes = data.ndata;
  uint nfeats = data.nfeats;
  const arr &psi = data.psi;
  const arr &psiT = data.psiT;
  arr &a = data.a;
  arr &b = data.b;
  double &sumloga = data.sumloga;
  double &sumlogb = data.sumlogb;
  arr &p = data.p;
  double &logz = data.logz;

  a.resize(nframes, 2);
  b.resize(nframes, 2);
  p.resize(nframes, 2);

  // in CRFs, a and b represent the message passing algorithm, not probability distributions!
  for(uint i = 0; i < 2; i++)
    a(0, i) = psi(0, i, 0);
  b[nframes - 1] = 1;
  sumloga = log(normalizeDist(a[0]()));
  sumlogb = log(normalizeDist(b[nframes - 1]()));
  for(uint ta = 1, tb = nframes - 1; ta < nframes; ta++, tb--) {
    a[ta]() = psi[ta] * a[ta-1];
    b[tb-1]() = psiT[tb] * b[tb];
    sumloga += log(normalizeDist(a[ta]()));
    sumlogb += log(normalizeDist(b[tb-1]()));
  }

  // CHECK(fabs(sumloga - sumlogb) < 1e-2, "sumloga (" << sumloga << ") and sumlogb (" << sumlogb << ") should be equal.");
  logz = .5 * (sumloga + sumlogb);

  p = a % b;
  for(uint t = 0; t < nframes; t++)
    normalizeDist(p[t]());
}

void LinChainCRF::update_pp(DataCRF &data, const arr &_theta) {
  update_p(data, _theta);

  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");
  if(data.last_theta_pp.N != 0 && data.last_theta_pp == theta)
    return;
  data.last_theta_pp = theta;

  uint nframes = data.ndata;
  uint nfeats = data.nfeats;
  const arr &a = data.a;
  const arr &b = data.b;
  const arr &phi = data.phi;
  const arr &psi = data.psi;
  const arr &p = data.p;
  arr &pp = data.pp;
  arr &dlogz = data.dlogz;

  pp.resize(nframes, 2, 2);
  dlogz.resize(nfeats);
  dlogz.setZero();

  for(uint i = 0; i < 2; i++)
    pp[0][i]() = p(0, i);
  for(uint t = 1; t < nframes; t++) {
    pp[t]() = (b[t] ^ a[t-1]) % psi[t];
    // for(uint i = 0; i < 2; i++)
    //   for(uint j = 0; j < 2; j++)
    //     pp(t, i, j) = a(t-1, j) * psi(t, i, j) * b(t, i);
    normalizeDist(pp[t]());
  }

  // for(uint t = 0; t < 10; t++)
  //   cout << "pp[" << t << "]: " << pp[t] << endl;

  for(uint t = 0; t < nframes; t++)
    for(uint i = 0; i < 2; i++)
      for(uint j = 0; j < 2; j++)
        dlogz += pp(t, i, j) * phi[t][i][j];
}

void LinChainCRF::update_ppp(DataCRF &data, const arr &_theta) {
  update_pp(data, _theta);

  const arr &theta = (&_theta != nullptr)? _theta: *params.get<arr>("theta");
  if(data.last_theta_ppp.N != 0 && data.last_theta_ppp == theta)
    return;
  data.last_theta_ppp = theta;

  uint nframes = data.ndata;
  const arr &psi = data.psi;
  const arr &a = data.a;
  const arr &b = data.b;
  const arr &pp = data.pp;
  arr &ppp = data.ppp;

  ppp.resize(TUP(nframes, 2, 2, 2));

  for(uint i = 0; i < 2; i++) {
    ppp[0][i]() = pp(0, i, 0);
    for(uint j = 0; j < 2; j++)
      ppp[1][i][j]() = pp(1, i, j);
  }
  for(uint t = 2; t < nframes; t++) {
    for(uint i = 0; i < 2; i++)
      for(uint j = 0; j < 2; j++)
        for(uint k = 0; k < 2; k++)
          ppp(TUP(t, i, j, k)) = a(t-2, k) * psi(t-1, j, k) * psi(t, i, j) * b(t, i);
    normalizeDist(ppp[t]());
  }
}

Model::Data *LinChainCRF::newData() {
  return new DataCRF();
}

void LinChainCRF::setFeatures(Data &data, const arr &d) {
  uint nframes = data.ndata;
  uint nfeats = data.nfeats;

  String transitions = MT::getParameter<String>("transitions");
  if(transitions == "true")
    nfeats += 4;
  data.phi.resize(TUP(nframes, 2, 2, nfeats));
  data.phi.setZero();
  data.nfeats = nfeats;
  params.set<uint>("nfeats", nfeats);

  if(transitions == "true") {
    for(uint t = 0; t < nframes; t++)
      for(uint i = 0; i < 2; i++)
        for(uint j = 0; j < 2; j++) {
          if(t != 0)
            data.phi[t][i][j](2 * i + j) = 1;
          if(i == 1)
            data.phi[t][i][j].subRange(4, -1) = d[t];
        }
  }
  else if(transitions == "false")
    for(uint t = 0; t < nframes; t++)
      for(uint i = 0; i < 2; i++)
        for(uint j = 0; j < 2; j++)
          if(i == 1)
            data.phi[t][i][j]() = d[t];
}

double LinChainCRF::negloglikelihood(const arr &x, Data &_data, const arr &_theta) {
  return 0./0.;
}

double LinChainCRF::negloglikelihood(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_pp(data, _theta);

  uint nframes = data.ndata;
  const arr &target = data.target;
  const arr &phiTtheta = data.phiTtheta;
  double &logz = data.logz;

  double nll;
  nll = -phiTtheta(0, target(0), 0);
  for(uint t = 1; t < nframes; t++)
    nll -= phiTtheta(t, target(t), target(t-1));
  nll += logz;

  return nll;
}

arr LinChainCRF::jacobian(Data &_data, const arr &_theta) { DataCRF &data = static_cast<DataCRF&>(_data);
  update_pp(data, _theta);

  uint nframes = data.ndata;
  uint nfeats = data.nfeats;
  const arr &target = data.target;
  const arr &phi = data.phi;
  const arr &dlogz = data.dlogz;

  arr J(nfeats);
  J = -phi[0][target(0)][0];
  for(uint t = 1; t < nframes; t++)
    J -= phi[t][target(t)][target(t-1)];
  J += dlogz;

  return J;
}

arr LinChainCRF::hessian(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_ppp(data, _theta);

  uint nframes = data.ndata;
  uint nfeats = data.nfeats;
  const arr &phi = data.phi;
  const arr &pp = data.pp;
  const arr &ppp = data.ppp;

  arr H(nfeats, nfeats);
  H.setZero();
  for(uint t = 0; t < nframes; t++)
    for(uint i = 0; i < 2; i++)
      for(uint j = 0; j < 2; j++)
        H += pp(t, i, j) * (phi[t][i][j] ^ phi[t][i][j]);
  for(uint t = 1; t < nframes; t++)
    for(uint i = 0; i < 2; i++)
      for(uint j = 0; j < 2; j++)
        for(uint k = 0; k < 2; k++)
          H += ppp(TUP(t, i, j, k)) * ((phi[t][i][j] ^ phi[t-1][j][k]) + (phi[t-1][j][k] ^ phi[t][i][j]));
  for(uint t = 0; t < nframes; t++)
    for(uint i = 0; i < 2; i++)
      for(uint j = 0; j < 2; j++)
        for(uint k = 0; k < 2; k++)
          for(uint l = 0; l < 2; l++)
            H -= pp(t, i, j) * pp(t, k, l) * (phi[t][i][j] ^ phi[t][k][l]);
  for(uint t = 1; t < nframes; t++)
    for(uint i = 0; i < 2; i++)
      for(uint j = 0; j < 2; j++)
        for(uint k = 0; k < 2; k++)
          for(uint l = 0; l < 2; l++)
            H -= pp(t, i, j) * pp(t-1, k, l)
                * ((phi[t][i][j] ^ phi[t-1][k][l]) + (phi[t-1][k][l] ^ phi[t][i][j]));
  return H;
}

arr LinChainCRF::sample(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_pp(data, _theta);

  double nframes = data.ndata;
  const arr &p = data.p;
  const arr &pp = data.pp;

  arr sample(nframes);
  sample(0) = (rnd.uni() <= p(0, 1));
  // cout << "p: " << p(0, 1) << " -> " << sample(0) << endl;
  for(uint t = 1; t < nframes; t++) {
    double condp = pp(t, 1, sample(t-1)) / p(t-1, sample(t-1));
    // sample(t) = (rnd.uni() <= pp(t, 1, sample(t-1))/p(t-1, sample(t-1)));
    sample(t) = (rnd.uni() <= condp);
    // cout << "sample(t-1): " << sample(t-1) << endl;
    // cout << "pp: " << pp(t, 1, sample(t-1)) << endl;
    // cout << "p: " << p(t-1, sample(t-1)) << endl;
    // cout << "condp: " << condp << endl;
    // cout << "sample(t): " << sample(t) << endl;
    // if(condp > 1)
    //   exit(1);
  }
  // exit(1);
  return sample;
}

arr LinChainCRF::maxlikelihood(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_psi(data, _theta);

  double nframes = data.ndata;
  const arr &phiTtheta = data.phiTtheta;

  arr vtmp(2);
  arr v(nframes + 1, 2);
  uintA vidx(nframes, 2);
  v.setZero();
  for(int t = nframes - 1; t >= 0; t--) {
    for(uint j = 0; j < 2; j++) {
      for(uint i = 0; i < 2; i++)
        vtmp(i) = phiTtheta(t, i, j) + v(t+1, i);
      vidx(t, j) = vtmp.maxIndex();
      v(t, j) = vtmp.max();
    }
  }

  arr vit(nframes);
  vit(0) = vidx(0, 0);
  for(uint t = 1; t < nframes; t++)
    vit(t) = vidx(t, vit(t-1));
  // cout << "vit: " << vit << endl;
  return vit;
}

arr LinChainCRF::marginal(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_p(data, _theta);
  return data.p;
}

arr LinChainCRF::marginal2(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_pp(data, _theta);
  return data.pp;
}

double LinChainCRF::nseq_expected(Data &_data, const arr &_theta) {
  DataCRF &data = static_cast<DataCRF&>(_data);
  update_pp(data, _theta);

  double nframes = data.ndata;
  const arr &pp = data.pp;

  double nseq = pp(0, 1, 0);
  for(uint t = 1; t < nframes; t++)
    nseq += pp(t, 1, 0);

  return nseq;
}

