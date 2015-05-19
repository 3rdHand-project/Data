#include <Optim/optimization.h>
#include <Gui/opengl.h>
#include <Gui/plot.h>
#include <Mocap/mocapgui.h>
#include <Mocap/utilAB.h>
#include "kf.h"
// #include "gui.h"
#include "optim.h"
#include "crf.h"
#include "logreg.h"

#define VOTE_MODE_STR(elem) case elem: return #elem;
const char *vote_mode_str(VoteMode mode) {
  switch(mode) {
    VOTE_MODE_LIST(VOTE_MODE_STR)
    default:
      HALT("VoteMode not defined?");
  }
};

KeyValueGraph votemode_kvg;
#define VOTE_MODE_TYPE(elem) votemode_kvg.append(#elem, new VoteMode(elem));
VoteMode str_to_votemode(const char *str) {
  if(votemode_kvg.N == 0) {
    VOTE_MODE_LIST(VOTE_MODE_TYPE);
  }
  VoteMode *mode = votemode_kvg.getValue<VoteMode>(str);
  return mode == nullptr? NO_VOTE_MODE: *mode;
}

// TODO reorder this
setf_t polyset(setf_t monoset) {
  return [&monoset] (const arr &inter) {
    uintA kflist;
    if(inter.nd == 1)
      kflist.append(monoset(inter));
    else if(inter.nd == 2)
      for(uint i = 0; i < inter.d0; i++)
        kflist.append(monoset(inter[i]));
    else if(inter.nd == 3)
      for(uint i = 0; i < inter.d0; i++)
        for(uint j = 0; j < inter.d1; j++)
          kflist.append(monoset(inter[i][j]));
    else HALT("Case not implemented yet.");
    std::sort(kflist.begin(), kflist.end());
    return kflist;
  };
}


// extern MocapWhitener g4w;

// KF_Score {{{
KF_Score::KF_Score() { clear(); }
KF_Score::KF_Score(const arr &t, const arr &y) {
  CHECK(t.N == y.N, "Output and target arrays must have the same size.");
  clear();
  tot = t.N;
  totpos = sum(t);
  totneg = tot - totpos;
  // TODO change types of target and y to uintA..
  for(uint i = 0; i < tot; i++)
    scores[(uint)y(i)][(uint)t(i)]++;
}

void KF_Score::clear() {
  tp = fp = fn = tp_ss = 0;
  tot = totpos = totneg = 0;
  memset(scores, 0, 4 * sizeof(uint));
}

double KF_Score::precision() const { return (double)tp/(double)(tp+fp); }
double KF_Score::recall() const { return (double)tp/(double)(tp+fn); }
double KF_Score::std() const { return sqrt((double)tp_ss/tp); }

// double KF_Score::tpr() const { return double(scores[1][1]) / totpos; }
// double KF_Score::fpr() const { return double(scores[1][0]) / totneg; }
// double KF_Score::accuracy() const { return double(scores[0][0] + scores[1][1]) / tot; }
// double KF_Score::f1() const { return 2.*scores[1][1] / (2.*scores[1][1] + scores[1][0] + scores[0][1]); }

KF_Score KF_Score::operator+(const KF_Score &x) {
  KF_Score s;
  s.tp = tp + x.tp;
  s.fp = fp + x.fp;
  s.fn = fn + x.fn;
  s.tp_ss = tp_ss + x.tp_ss;

  // s.tot = tot + x.tot;
  // s.totpos = totpos + x.totpos;
  // s.totneg = totneg + x.totneg;
  // s.scores[0][0] = scores[0][0] + x.scores[0][0];
  // s.scores[0][1] = scores[0][1] + x.scores[0][1];
  // s.scores[1][0] = scores[1][0] + x.scores[1][0];
  // s.scores[1][1] = scores[1][1] + x.scores[1][1];
  return s;
}

KF_Score& KF_Score::operator+=(const KF_Score &x) {
  tp += x.tp;
  fp += x.fp;
  fn += x.fn;
  tp_ss += x.tp_ss;
  // tot += x.tot;
  // totpos += x.totpos;
  // totneg += x.totneg;
  // scores[0][0] += x.scores[0][0];
  // scores[0][1] += x.scores[0][1];
  // scores[1][0] += x.scores[1][0];
  // scores[1][1] += x.scores[1][1];
  return *this;
}

void KF_Score::write(std::ostream &os) const {
  os << " tp / fp / fn: " << tp << " / " << fp << " / " << fn << endl;
  os << " precision / recall: " << precision() << " / " << recall() << endl;
  os << " std: " << std() << endl;
  // os << "tot: " << tot << " (neg: " << totneg << ", pos: " << totpos << ")" << endl;
  // os << "     t=0    t=1" << endl;
  // os << "y=0 " << std::setw(5) << std::fixed << std::setprecision(2) << (1-fpr())
  //     << " " << std::setw(5) << std::fixed << std::setprecision(2) << (1-tpr()) << endl;
  // os << "y=1 " << std::setw(5) << std::fixed << std::setprecision(2) << fpr()
  //     << " " << std::setw(5) << std::fixed << std::setprecision(2) << tpr() << endl;
  // os << "accuracy: " << accuracy() << endl;
  // os << "f1: " << f1() << endl;
}
// KF_Score }}}
// ROCCurve {{{
// void ROCCurve::append(double thresh, const KF_Score &score) {
//   thresholds.append(thresh);
//   scores.append(score);
//   ind.append(0);
// }
// double ROCCurve::area() {
//   sort();

//   double tpr1, tpr2, fpr1, fpr2;
//   double a = 0;
//   for(uint i = 0; i < ind.N -1; i++) {
//     tpr1 = scores(ind(i)).tpr();
//     tpr2 = scores(ind(i+1)).tpr();
//     fpr1 = scores(ind(i)).fpr();
//     fpr2 = scores(ind(i+1)).fpr();
//     a += .5 * (fpr2 - fpr1) * (tpr1 + tpr2);
//   }
//   a += .5 * scores(ind(0)).fpr() * scores(ind(0)).tpr()
//        + .5 * (1 - scores(ind.last()).fpr()) * (1 + scores(ind.last()).tpr());
//   return a;
// }

// void ROCCurve::show(const char *pdf) {
//   sort();

//   FGP fgplot, fgplot_pdf;
//   fgplot.set("domain", true);
//   fgplot.set("lines", true);
//   fgplot.set("ymin_b", true);
//   fgplot.set("ymax_b", true);
//   fgplot.set("ymin", 0.);
//   fgplot.set("ymax", 1.);
//   fgplot.open();
//   if(pdf) {
//     fgplot_pdf.set("domain", true);
//     fgplot_pdf.set("lines", true);
//     fgplot_pdf.set("ymin_b", true);
//     fgplot_pdf.set("ymax_b", true);
//     fgplot_pdf.set("ymin", 0.);
//     fgplot_pdf.set("ymax", 1.);
//     fgplot_pdf.set("hardcopy", pdf);
//     fgplot_pdf.open();
//   }
//   for(uint i = 0; i < scores.N; i++) {
//     fgplot() << scores(ind(i)).fpr() << " " << scores(ind(i)).tpr();
//     if(pdf)
//       fgplot_pdf() << scores(ind(i)).fpr() << " " << scores(ind(i)).tpr();
//   }
//   fgplot.close();
//   if(pdf)
//     fgplot_pdf.close();
// }
// void ROCCurve::sort() {
//   if(scores.N > 1 && ind.last() == 0) {
//     for(uint i = 0; i < ind.N; i++)
//       ind(i) = i;
//     std::sort(ind.begin(), ind.end(),
//               [this](uint i1, uint i2) {
//                 return scores(i1).fpr() < scores(i2).fpr();
//               });
//   }
// }
// }}}
// KF_Model {{{
// KF_Model::KF_Model(Target target): mod(nullptr), target(target) { target_str = Target_to_str(target); setDefaultParams(); }
KF_Model::KF_Model(Target target): mod(nullptr), target(target), target_str(Target_to_str(target)) { setDefaultParams(); }
KF_Model::~KF_Model() { if(mod) delete mod; }

Model &KF_Model::model() {
  if(mod == nullptr) {
    String model = *params.get<String>("model");
    if(model == "crf")
      mod = new LinChainCRF();
    else if(model == "logreg")
      mod = new LogReg();
    else HALT("No model called '" << model << "'.");
  }
  return *mod;
}

void KF_Model::setData(Model::Data &data, MocapSeq &seq) {
  data.ndata = seq.nframes_thin;
  data.nfeats = seq.data.d1;
  data.target = seq.ann_thin;

  model().setFeatures(data, seq.data);
}

void KF_Model::setDefaultParams() {
  String model = MT::getParameter<String>("model");
  uint nsplits = MT::getParameter<double>("nsplits");
  String output = MT::getParameter<String>("output");
  String vote_mode = MT::getParameter<String>("vote_mode");

  params.set("model", model);
  params.set("nsplits", nsplits);
  params.set("output", output);
  params.set("vote_mode", vote_mode);
};

arr KF_Model::thinnen(double nframes_thin, const arr &thick) {
  arr thin(nframes_thin);
  uint thinning = MT::getParameter<uint>("thinning");
  for(uint f_thin = 0; f_thin < nframes_thin; f_thin++)
    thin(f_thin) = thick((f_thin+1) * thinning - 1);

  return thin;
}

arr KF_Model::thicken(double nframes, const arr &thin) {
  uint nframes_thin = thin.N;
  arr thick(nframes);
  uint thinning = MT::getParameter<uint>("thinning");
  for(uint f_thin = 0; f_thin < nframes_thin; f_thin++)
    thick.subRange(f_thin * thinning, (f_thin + 1) * thinning - 1)() = thin(f_thin);
  if(nframes_thin * thinning < nframes)
    thick.subRange(nframes_thin * thinning, -1)() = thin.last();

  return thick;
}

arr KF_Model::maxlikelihood(MocapRec &mrec, const String &limb, const String &part, Thickness thickness) {
  const StringA &digits = mrec.id().sensorsof(limb);
  uint ndigits = digits.N;

  // TODO set number of frames....
  uint nframes_thin = mrec.numFrames(THIN);
  arr interactions(ndigits, nframes_thin);
  for(uint di = 0; di < ndigits; di++) {
    MocapSeq *seq = mrec.seq(digits(di), part);
    interactions[di]() = query(*seq, THIN);
    delete seq;
  }

  arr interaction = sum(interactions, 0);

  String vote_mode = *params.get<String>("vote_mode");
  if(vote_mode == "majority") {
    for(uint i = 0; i < interaction.N; i++)
      interaction(i) = (interaction(i) >= ndigits/2.)? 1.: 0.;
  }
  else if(vote_mode == "unanimous") {
    // should evaluate to 1 only if full vote.
    for(uint i = 0; i < interaction.N; i++)
      interaction(i) = (interaction(i) == ndigits)? 1.: 0.;
  }
  else if(vote_mode == "nega-unanimous") {
    for(uint i = 0; i < interaction.N; i++)
      interaction(i) = (interaction(i) > 0)? 1.: 0.;
  }
  else HALT("Vote_mode ''" << vote_mode << "' not programmed.");

  if(thickness == THICK)
    return thicken(mrec.numFrames(THICK), interaction);
  return interaction;
}

arr KF_Model::marginal(MocapRec &mrec, const String &limb, const String &part) {
  CHECK(mrec.id().agent_limbs().contains(limb), "Recording does not contain limb '" << limb << "'");
  CHECK(mrec.id().object_parts().contains(part), "Recording does not contain part '" << part << "'");

  uint nframes_thin = mrec.numFrames(THIN);
  const StringA &digits = mrec.id().sensorsof(limb);
  uint ndigits = digits.N;

  // arr marginals(ndigits, nframes_thin, 2);
  uint nsamples = 1000;
  arr cumsamples(nsamples, nframes_thin);
  for(uint di = 0; di < ndigits; di++) {
    MocapSeq &seq = *mrec.seq(digits(di), part);
    Model::Data &data = *model().newData();
    setFeats(seq);
    setData(data, seq);
    // marginals[di]() = model().marginal(data);
    for(uint si = 0; si < nsamples; si++)
      cumsamples[si]() += model().sample(data);
    delete &seq;
  }

  // TODO Generalize to generic dimensions
  // TODO Improve sum......
  // for(uint t = 0; t < nframes_thin; t++) {
  //   for(uint i = 0; i < 2; i++)
  //     for(uint j = 0; j < 2; j++)
  //       for(uint k = 0; k < 2; k++)
  //         if(i+j+k>1)
  //           marginal(t, 1) += marginals(0, t, i)
  //                           * marginals(1, t, j)
  //                           * marginals(2, t, k);
  //   marginal(t, 0) = 1-marginal(t, 1);
  // }
  // for(uint i = 0; i < 100; i++)
  //   marginal += model().sample()

  arr marginal(nframes_thin, 2);
  marginal.setZero();
  for(uint t = 0; t < nframes_thin; t++)
    for(uint si = 0; si < nsamples; si++)
      if(cumsamples(si, t) > 1)
        marginal(t, 1)++;
      else
        marginal(t, 0)++;
  for(uint t = 0; t < nframes_thin; t++)
    normalizeDist(marginal[t]());

  return marginal;
}

arr KF_Model::marginal2(MocapRec &mrec, const String &limb, const String &part) {
  CHECK(mrec.id().agent_limbs().contains(limb), "Recording does not contain limb '" << limb << "'");
  CHECK(mrec.id().object_parts().contains(part), "Recording does not contain part '" << part << "'");

  uint nframes_thin = mrec.numFrames(THIN);
  const StringA &digits = mrec.id().sensorsof(limb);
  uint ndigits = digits.N;

  arr marginal(nframes_thin, 2, 2);
  arr marginals; marginals.resize(TUP(ndigits, nframes_thin, 2, 2));
  for(uint di = 0; di < ndigits; di++) {
    MocapSeq &seq = *mrec.seq(digits(di), part);
    Model::Data &data = *model().newData();
    setFeats(seq);
    setData(data, seq);
    marginals[di]() = model().marginal2(data);
    delete &seq;
  }

  // TODO Generalize to generic dimensions
  // TODO Improve sum......
  marginal.setZero();
  for(uint t = 0; t < nframes_thin; t++)
    for(uint l = 0; l < 2; l++)
    for(uint ll = 0; ll < 2; ll++)
    for(uint i = 0; i < 2; i++)
    for(uint ii = 0; ii < 2; ii++)
      for(uint j = 0; j < 2; j++)
      for(uint jj = 0; jj < 2; jj++)
        for(uint k = 0; k < 2; k++)
        for(uint kk = 0; kk < 2; kk++)
          if(((l && i+j+k>1) || (!l && i+j+k<2))
            && ((ll && ii+jj+kk>1) || (!ll && ii+jj+kk<2)))
            marginal(t, l, ll) += marginals(TUP(0, t, i, ii))
                                * marginals(TUP(1, t, j, jj))
                                * marginals(TUP(2, t, k, kk));

  return marginal;
}

arr KF_Model::entropy(MocapRec &mrec, const String &limb, const String &part) {
  arr marg = marginal(mrec, limb, part);
  uint nframes = marg.d0;

  arr H(nframes);
  for(uint t = 0; t < marg.d0; t++) {
    H(t) = -sum(marg[t]%log(marg[t]));
    if(std::isnan(H(t))) H(t) = 0;
  }
  return H;
}

uint KF_Model::nseq_ground_truth(MocapRec &mrec, const String &limb, const String &part) {
  // arr gt = mrec.ann(target, mrec.id().sensorsof(limb).elem(0), part);
  const String &digit = mrec.id().sensorsof(limb).elem(0);
  arr &gt = *mrec.label().getValue<arr>(STRINGS(target_str, digit, part));
  if(!gt.N) return 0;

  uint nseq = gt(0);
  for(uint t = 1; t < gt.d0; t++)
    if(gt(t) && !gt(t-1))
      nseq++;
  return nseq;
}

uint KF_Model::nseq_maxlikelihood(MocapRec &mrec, const String &limb, const String &part) {
  arr ml = maxlikelihood(mrec, limb, part, THIN);
  uint nframes_thin = mrec.numFrames(THIN);

  uint nseq = ml(0);
  for(uint t = 1; t < nframes_thin; t++)
    if(ml(t) && !ml(t-1))
      nseq++;
  return nseq;
}

double KF_Model::nseq_expected(MocapRec &mrec, const String &limb, const String &part) {
  arr pp = marginal2(mrec, limb, part);
  uint nframes_thin = mrec.numFrames(THIN);

  double nseq = pp(0, 1, 0);
  for(uint t = 1; t < nframes_thin; t++)
    nseq += pp(t, 1, 0);

  return nseq;
}

uintA KF_Model::lseq_ground_truth(MocapRec &mrec, const String &limb, const String &part) {
  const String &digit = mrec.id().sensorsof(limb).elem(0);
  arr &gt = *mrec.label().getValue<arr>(STRINGS(target_str, digit, part));
  // arr gt = mrec.ann(target, mrec.id().sensorsof(limb).elem(0), part);
  uintA lseq;

  if(!gt.N) return lseq;

  if(gt(0))
    lseq.append(1);
  for(uint t = 1; t < gt.d0; t++)
    if(gt(t)) {
      if(!gt(t-1))
        lseq.append(1);
      else
        lseq.last()++;
    }
  return lseq;
}

uintA KF_Model::lseq_maxlikelihood(MocapRec &mrec, const String &limb, const String &part) {
  arr ml = maxlikelihood(mrec, limb, part, THICK);
  uint nframes = mrec.numFrames(THICK);

  uintA lseq;
  if(ml(0))
    lseq.append(1);
  for(uint t = 1; t < nframes; t++)
    if(ml(t)) {
      if(!ml(t-1))
        lseq.append(1);
      else
        lseq.last()++;
    }
  return lseq;
}

uintA KF_Model::lseq_expected(MocapRec &mrec, const String &limb, const String &part) {
  return uintA();
}

arr KF_Model::sample(MocapRec &mrec, const String &limb, const String &part) {
  CHECK(mrec.id().agent_limbs().contains(limb), "Recording does not contain limb '" << limb << "'");
  CHECK(mrec.id().object_parts().contains(part), "Recording does not contain part '" << part << "'");

  uint nframes_thin = mrec.numFrames(THIN);
  arr sample(nframes_thin);
  sample.setZero();

  const StringA &digits = mrec.id().sensorsof(limb);
  uint ndigits = digits.N;
  for(uint di = 0; di < ndigits; di++) {
    MocapSeq &seq = *mrec.seq(digits(di), part);
    Model::Data &data = *model().newData();

    setFeats(seq);
    setData(data, seq);

    sample += model().sample(data);

    delete &seq;
    delete &data;
  }
  for(uint t = 0; t < nframes_thin; t++)
    sample(t) = (sample(t)>1? 1: 0);
  return sample;
}

double KF_Model::precision(MocapRec &mrec, const String &limb, const String &part) {
  arr ml = maxlikelihood(mrec, limb, part, THICK);
  // arr gtruth = mrec.ann(target, limb, part);
  arr &gtruth = *mrec.label().getValue<arr>(STRINGS(target_str, limb, part));

  double precision, total;
  precision = total = 0;
  for(uint t = 0; t < gtruth.d0; t++) {
    if(ml(t)) {
      total++;
      if(gtruth(t))
        precision++;
    }
  }
  return double(precision) / double(total);
}

double KF_Model::recall(MocapRec &mrec, const String &limb, const String &part) {
  arr ml = maxlikelihood(mrec, limb, part, THICK);
  // arr gtruth = mrec.ann(target, limb, part);
  arr &gtruth = *mrec.label().getValue<arr>(STRINGS(target_str, limb, part));

  double recall, total;
  recall = total = 0;
  for(uint t = 0; t < gtruth.d0; t++) {
    if(gtruth(t)) {
      total++;
      if(ml(t))
        recall++;
    }
  }
  return (double)recall / (double)total;
}

double KF_Model::estimate_entropy(MocapRec &mrec, const String &limb, const String &part) {
  // CHECK(mrec.id().agent_limbs().contains(limb), "Recording does not contain limb '" << limb << "'");
  // CHECK(mrec.id().object_parts().contains(part), "Recording does not contain part '" << part << "'");

  // double entropy = 0;

  // uint nframes_thin = mrec.numFrames(THIN);
  // arr sample(nframes_thin);
  // sample.setZero();

  // const StringA &digits = mrec.id().sensorsof(limb);
  // uint ndigits = digits.N;

  // for(uint di = 0; di < ndigits; di++) {
  //   MocapSeq &seq = *mrec.seq(digits(di), part);
  //   Model::Data &data = *model().newData();

  //   setFeats(seq);
  //   setData(data, seq);

  //   s = model().sample(data);
  //   nll = model().neglikelihood(data);
  //   sample += s;

  //   delete &seq;
  //   delete &data;
  // }
  // return entropy;
  return 0./0.;
}

void KF_Model::plot_profile(MocapRec &mrec, const String &limb, const String &part) {
  CHECK(mrec.id().agent_limbs().contains(limb), "Recording does not contain limb '" << limb << "'");
  CHECK(mrec.id().object_parts().contains(part), "Recording does not contain part '" << part << "'");

  const StringA &digits = mrec.id().sensorsof(limb);
  uint ndigits = digits.N;

  plotClear();
  for(uint di = 0; di < ndigits; di++) {
    MocapSeq &seq = *mrec.seq(digits(di), part);
    Model::Data &data = *model().newData();

    setFeats(seq);
    setData(data, seq);

    plotFunction(model().marginal(data).col(1));

    delete &seq;
    delete &data;
  }
  arr marg = marginal(mrec, limb, part).col(1);
  arr H = entropy(mrec, limb, part);
  arr ml = maxlikelihood(mrec, limb, part, THIN);
  arr ml_thick = maxlikelihood(mrec, limb, part, THICK);
  // arr ann = mrec.ann(target, limb, part);
  arr &ann = *mrec.label().getValue<arr>(STRINGS(target_str, limb, part));
  arr gtruth = thinnen(mrec.numFrames(THIN), ann);

  if(mrec.dir.contains(STRING("coop_1")) && part == "/toolbox/side_left") {
    // arr tmp;
    // tmp.append(ann);
    // tmp.append(ml_thick);
    // tmp.reshape(2, tmp.N/2);
    // tmp = ~tmp;

    // cout << "=== " << limb << " - " << part << " ===" << endl;
    // uint last = 0;
    // cout << tmp[0] << " @ " << 0 << endl;
    // for(uint i = 0; i < tmp.d0; i++) {
    //   if(tmp[i] != tmp[last]) {
    //     cout << tmp[i] << " @ " << i << " / " << tmp.d0 << endl;
    //     last = i;
    //   }
    // }
    cout << "=== " << limb << " - " << part << " ===" << endl;
    cout << "ann: ";
    if(ann(0))
      cout << "0 - ";
    for(uint t = 1; t < ann.N; t++) {
      if(ann(t) && !ann(t-1))
        cout << (t-2100.)/550. << " - ";
      if(!ann(t) && ann(t-1))
        cout << (t-1-2100.)/550. << ", ";
    }
    cout << endl;

    cout << "ml_thick: ";
    if(ml_thick(0))
      cout << "0 - ";
    for(uint t = 1; t < ml_thick.N; t++) {
      if(ml_thick(t) && !ml_thick(t-1))
        cout << (t-2100.)/550. << " - ";
      if(!ml_thick(t) && ml_thick(t-1))
        cout << (t-1-2100.)/550. << ", ";
    }
    cout << endl;
  }

  // plotFunction(marg);
  plotFunction(H);
  plotFunction(ml);
  plotFunction(gtruth);
  // cout << "marg.getDim(): " << marg.getDim() << endl;
  // cout << "ml.getDim(): " << ml.getDim() << endl;
  // cout << "ml: " << ml << endl;

  // plotFunction(marginal(mrec, limb, part).col(1));
  // plotFunction(maxlikelihood(mrec, limb, part));
  plot();
}

void KF_Model::plot_samples(MocapRec &mrec, const String &limb, const String &part) {
  CHECK(mrec.id().agent_limbs().contains(limb), "Recording does not contain limb '" << limb << "'");
  CHECK(mrec.id().object_parts().contains(part), "Recording does not contain part '" << part << "'");

  const StringA &digits = mrec.id().sensorsof(limb);
  uint ndigits = digits.N;
  for(uint di = 0; di < ndigits; di++) {
    MocapSeq &seq = *mrec.seq(digits(di), part);
    Model::Data &data = *model().newData();

    setFeats(seq);
    setData(data, seq);

    for(uint i = 0; i < 10; i++)
      plotFunction(model().sample(data));
    plot();
    exit(1);

    delete &seq;
    delete &data;
  }
  plot();
}

KF_Score KF_Model::statistics(MocapRec &mrec) {
  KF_Score score;

  StringA elist1, elist2;
  entitylists(mrec, elist1, elist2);
  for(const String &entity1: elist1)
    for(const String &entity2: elist2)
      score += statistics(mrec, entity1, entity2);

  return score;
}

KF_Score KF_Model::statistics(MocapRec &mrec, const String &o1, const String &o2) {
  KF_Score score;

//   mrec.
//   gt_label = mrec.ann(STRINGS(Target_to_str(target), o1, o2));
//   ml_label = maxlikelihood(mrec, o1, o2, THICK);

  arr gt, ml;
  uint nframes = mrec.numFrames(THICK) + 2;

  // TODO check annotation for target, o1, o2 pair
  gt.append(0);
  // gt.append(mrec.ann(STRINGS(Target_to_str(target), o1, o2)));
  arr &tmpgt = *mrec.label().getValue<arr>(STRINGS(target_str, o1, o2));
  gt.append(tmpgt);
  gt.append(0);

  ml.append(0);
  ml.append(maxlikelihood(mrec, o1, o2, THICK));
  ml.append(0);

  // TODO
  // KF_Score statistics(const arr &seq1, const arr &seq2);

  // CHECK(seq1.nd == 1, "First sequence has too many dimensions.");
  // CHECK(seq2.nd == 1, "Second sequence has too many dimensions.");
  // CHECK(seq1.N == seq2.N, "Sequence lengths don't match.");
  // KF_Score score;
  // arr gt, ml;
  // uint nframes = seq1.N;

  // gt.append(0);
  // gt.append(seq1);
  // gt.append(0);

  // ml.append(0);
  // ml.append(seq2);
  // ml.append(0);
  // nframes += 2;

  int gt_kf, ml_kf;
  uint ngt_kf, nml_kf;
  bool correct, correct_prev;

  gt_kf = ml_kf = -1;
  ngt_kf = nml_kf = 0;
  correct = true;
  for(uint t = 1; t < nframes; t++) {
    correct_prev = correct;
    correct = (ml(t) == gt(t));

    if(gt(t) != gt(t-1)) {
      ngt_kf++;
      gt_kf = (gt_kf == -1)? t: -1;
    }
    if(ml(t) != ml(t-1)) {
      nml_kf++;
      ml_kf = (ml_kf == -1)? t: -1;
    }
    if(correct && !correct_prev && gt_kf != -1 && ml_kf != -1) {
      score.tp++;
      score.tp_ss += (gt_kf - ml_kf) * (gt_kf - ml_kf);
      gt_kf = -1;
      ml_kf = -1;
    }
  }
  score.fp = nml_kf - score.tp;
  score.fn = ngt_kf - score.tp;
  // score.tpstd = sqrt(score.tpstd / score.tp);

  return score;
}

void KF_Model::print_statistics(MocapRec &mrec) {
  cout << "#########################" << endl;
  cout << "### " << mrec.dir << " ###" << endl;
  cout << "#########################" << endl << endl;
  cout << "### Digit statistics ###" << endl;
  for(const String &digit: mrec.id().agent_digits()) {
    for(const String &part: mrec.id().object_parts()) {
      cout << "# " << digit << " - " << part << endl;
      // arr ann = mrec.ann(Target::TS, digit, part);
      arr &ann = *mrec.label().getValue<arr>(STRINGS(Target_to_str(Target::TS), digit, part));
      if(ann.N) {
        double real_number = ann(0);
        for(uint t = 1; t < ann.N; t++)
          if(ann(t) && !ann(t-1))
            real_number++;

        MocapSeq &seq = *mrec.seq(digit, part);
        Model::Data &data = *model().newData();
        setFeats(seq);
        setData(data, seq);

        // TODO get sequence relative to digit-part pair
        cout << " - nseq_ground_truth: " << real_number << endl;
        cout << " - nseq_expected:     " << model().nseq_expected(data) << endl;
        cout << " - nseq_ml:           " << model().nseq_maxlikelihood(data) << endl;
        cout << " - entropy:           " << model().entropy(data) << endl;
      }
    }
  }
  cout << "### Limb statistics ###" << endl;
  for(const String &limb: mrec.id().agent_limbs()) {
    for(const String &part: mrec.id().object_parts()) {
      cout << "# " << limb << " - " << part << endl;
      cout << " - nseq_ground_truth:   " << nseq_ground_truth(mrec, limb, part) << endl;
      // cout << " - nseq_expected:       " << nseq_expected(mrec, limb, part) << endl;
      cout << " - nseq_ml:             " << nseq_maxlikelihood(mrec, limb, part) << endl;
      cout << " - lseq_ground_truth:   " << lseq_ground_truth(mrec, limb, part) << endl;
      cout << " - lseq_ml:             " << lseq_maxlikelihood(mrec, limb, part) << endl;
      cout << " - precision:           " << precision(mrec, limb, part) << endl;
      cout << " - recall:              " << recall(mrec, limb, part) << endl;
      cout << " - entropy:             " << estimate_entropy(mrec, limb, part) << endl;
    }
  }
}

// KF_Score KF_Model::test(const MocapSeqL &seqlist) {
//   KF_Score score;
//   for(MocapSeq *seq: seqlist)
//     score += KF_Score(seq->ann, query(*seq));
//   return score;
// }

// KF_Score KF_Model::crossval(const MocapSeqL &seqlist) {
//   KF_Score score;

//   uint nsplits = *params.get<uint>("nsplits");

//   MT::Array<MocapSeqL> trainseqlist, testseqlist;
//   split(trainseqlist, testseqlist, seqlist, nsplits);
//   for(uint n = 0; n < nsplits; n++) {
//     // take nth training split
//     train(trainseqlist(n));
//     // take nth testing split
//     score += test(testseqlist(n));
//   }
//   return score;
// }

// ROCCurve KF_Model::ROC(const MocapSeqL &seqlist, const arr &thresholds) {
//   ROCCurve roc;
//   KF_Score score;
//   ProgressBar pb;

//   pb.prefix << "Calc. ROC curve";
//   pb.reset(thresholds.N);
//   for(double thresh: thresholds) {
//     params.set("thresh", thresh);
//     score = crossval(seqlist);
//     pb.step() << "thresh = " << thresh;
//     roc.append(thresh, score);
//   }
//   pb.stop() << "done";
//   return roc;
// }

double KF_Model::negloglikelihood(MocapSeqL &seqlist) {
  double nll = 0;
  for(MocapSeq *seq: seqlist)
    nll += negloglikelihood(*seq);
  return nll;
}

// uintA extract_constraints(const StringA &objects) {
//   // TODO create pairs
//   StringA objectsensors;
//   for(const String &object: objects)
//     objectsensors.append(g4rec.id().sensorsof(object));
//   uint nagentsensors = objectsensors.N;
//   String &root = objectsensors(0);
// }
// }}}
// KF_Constraint {{{
KF_Constraint::KF_Constraint(const KF_Constraint &c) {
  main_obj = c.main_obj;
  obj = c.obj;
  frame = c.frame;
  p_mean = c.p_mean;
  q_mean = c.q_mean;
  p_cov = c.p_cov;
  q_cov = c.q_cov;
}

void KF_Constraint::write(std::ostream &os) const {
  os << "Constraint:" << endl;
  os << " * frame: " << frame << endl;
  os << " * main_obj: " << main_obj << endl;
  os << " * obj: " << obj << endl;
  os << " * p_mean: " << p_mean << endl;
  // os << " * p_cov: " << p_cov << endl;
  os << " * q_mean: " << q_mean << endl;
  // os << " * q_cov: " << q_cov << endl;
}
// }}}
// KF_CE_CRF {{{
// KFConstraintL KF_CE_CRF::extract(MocapRec &rec, KFModel &model, const StringA &objects) {
//   StringA objectsensors;
//   for(String &object: objects)
//     objectsensors.append(rec.id().sensorsof(object));

//   CHECK(objectsensors.N > 1, "Constraint detection requires at least 2 objects");

//   StringA tmpstrings;
//   arr tmpfs;
//   arr inter;
//   arr f(3);
//   String o1 = objectsensors(0);
//   String o2 = objectsensors(1);
//   inter = model.query(rec.seq(o1, o2));
//   f(0) = inter.findValue(1);
//   for(uint i = 2; i < objectsensors.N; i++) {
//     String &objectsensor = objectsensors(i);

//     inter = model.query(rec.seq(o1, objectsensor));
//     f(1) = inter.findValue(1);
//     inter = model.query(rec.seq(o2, objectsensor));
//     f(2) = inter.findValue(1);

//     switch(f.minIndex()) {
//       case 0:
//         tmpstrings.append(objectsensor);
//         tmpfs.append((f(1) + f(2)) / 2);
//         break;
//       case 1:
//         tmpstrings.append(o2);
//         tmpfs.append((f(0) + f(2)) / 2);

//         f(0) = f(1);
//         o2 = objectsensor;
//         break;
//       case 2:
//         tmpstrings.append(o1);
//         tmpfs.append((f(0) + f(1)) / 2);

//         f(0) = f(2);
//         o1 = objectsensor;
//         break;
//     }
//   }

//   KFConstraintL constraintl;
//   constraintl.append(computeConstraint(rec, o1, o2, f(0)));

//   return constraintl;
// }

// KFConstraint *KF_CE_CRF::computeConstraint(MocapRec &rec, const char *s1, const char *s2, uint f) {
//   KFConstraint *c = nullptr;
//   return c;
// }
// KFConstraintL KF_CE_CRF::findConstraints(MocapRec &rec, KFModel &model, const StringA &agents, const StringA &objects) {
//   const String &main_obj = objects(0);

//   arr inter, inter1, inter2;
//   MocapSeq seq;
//   uintA kflist;

//   KFConstraint *c;
//   KFConstraintL cl;
//   for(const String &obj: objects) {
//     if(obj == main_obj) continue;
//     cout << "===" << endl;
//     cout << main_obj << " - " << obj << endl;
//     seq = rec.seq(main_obj, obj);
//     inter = model.query(seq);
//     inter1 = model.maxlikelihood(rec, STRINGS("rh"), main_obj);
//     inter2 = model.maxlikelihood(rec, STRINGS("rh"), obj);

//     kflist.clear();
//     for(uint i = 1; i < inter.N; i++)
//       if(inter(i-1) != inter(i))
//         kflist.append(i);
//     for(uint i = 1; i < inter1.N; i++)
//       if(inter1(i-1) != inter1(i))
//       // if(inter1(i-1) && !inter1(i))
//         kflist.append(i);
//     for(uint i = 1; i < inter2.N; i++)
//       if(inter2(i-1) != inter2(i))
//       // if(inter2(i-1) && !inter2(i))
//         kflist.append(i);
//     std::sort(kflist.begin(), kflist.end());

//     c = new KFConstraint();
//     *c = findConstraint(rec, main_obj, obj, kflist);
//     cout << *c << endl;
//     cl.append(c);
//   }
//   return cl;
// }

// KFConstraint KF_CE_CRF::findConstraint(MocapRec &rec, const char *s1, const char *s2, const uintA &kflist) {
//   KFConstraint c, ct;
//   double l, lt;
  
//   // if(STRING(s2) == "toolbox:side_left"){
//   //   uint i;
//   //   for(i = 0; kflist(i) < 3100; i++);
//   //   c = computeConstraint(rec, s1, s2, kflist.subRange(i, -1));
//   //   c.frame = kflist(i);
//   //   c.main_obj = s1;
//   //   c.obj = s2;
//   // }
//   // if(STRING(s2) == "toolbox:side_back") {
//   //   uint i;
//   //   for(i = 0; kflist(i) < 7000; i++);
//   //   c = computeConstraint(rec, s1, s2, kflist.subRange(i, -1));
//   //   c.frame = kflist(i);
//   //   c.main_obj = s1;
//   //   c.obj = s2;
//   // }
//   // if(STRING(s2) == "toolbox:handle") {
//   //   uint i;
//   //   for(i = 0; kflist(i) < 10500; i++);
//   //   c = computeConstraint(rec, s1, s2, kflist.subRange(i, -1));
//   //   c.frame = kflist(i);
//   //   c.main_obj = s1;
//   //   c.obj = s2;
//   // }
//   // if(STRING(s2) == "toolbox:side_right") {
//   //   uint i;
//   //   for(i = 0; kflist(i) < 17000; i++);
//   //   c = computeConstraint(rec, s1, s2, kflist.subRange(i, -1));
//   //   c.frame = kflist(i);
//   //   c.main_obj = s1;
//   //   c.obj = s2;
//   // }
//   // return c;

//   bool first = true;
//   l = 0;
//   for(uint i = 0; i < kflist.N; i++) {
//     ct = computeConstraint(rec, s1, s2, kflist.subRange(i, -1));
//     lt = computeConstraintScore(rec, s1, s2, kflist, ct);
//     if(lt > l || first) {
//       first = false;
//       c = ct;
//       c.frame = kflist(i);
//       c.main_obj = s1;
//       c.obj = s2;
//       l = lt;
//     }
//   }
//   return c;
// }

// KFConstraint KF_CE_CRF::computeConstraint(MocapRec &rec, const char *s1, const char *s2, const uintA &kflist) {
//   KFConstraint c;
//   c.frame = kflist(0);

//   rec.computeDPos(s1);
//   rec.computeDQuat(s1);

//   arr dpos(kflist.N, 3), dquat(kflist.N, 4);
//   for(uint i = 0; i < kflist.N; i++) {
//     dpos[i]() = rec.query(STRING(s1 << "_dPos"), s2, kflist(i));
//     dquat[i]() = rec.query(STRING(s1 << "_dQuat"), s2, kflist(i));
//   }

//   c.p_mean = sum(dpos, 0) / (double)dpos.d0;
//   c.p_cov = zeros(3, 3);
//   for(uint i = 0; i < dpos.d0; i++)
//     for(uint j =0; j < dpos.d0; j++)
//       c.p_cov += dpos[i] ^ dpos[i];
//   c.p_cov /= (double)dpos.d0;
//   c.p_cov -= c.p_mean ^ c.p_mean;

//   c.q_mean = sum(dquat, 0) / (double)dquat.d0;
//   c.q_mean /= sqrt(sumOfSqr(c.q_mean));
//   c.q_cov = zeros(4, 4);
//   for(uint i = 0; i < dquat.d0; i++)
//     for(uint j =0; j < dquat.d0; j++)
//       c.q_cov += dquat[i] ^ dquat[i];
//   c.q_cov /= (double)dquat.d0;
//   c.q_cov -= c.q_mean ^ c.q_mean;

//   return c;
// }

// double KF_CE_CRF::computeConstraintScore(MocapRec &rec, const char *s1, const char *s2, const uintA &kflist, KFConstraint &c) {
//   double m = 2;

//   uint p_k = c.p_mean.N;
//   // c.p_cov = eye(3);
//   arr p_cov = c.p_cov;
//   arr p_covm = m * c.p_cov;
//   arr inv_p_cov = inverse(c.p_cov);
//   arr inv_p_covm = inverse(p_covm);
//   double p_det = determinant(c.p_cov);
//   double p_detm = determinant(p_covm);
//   double p_norm = 1. / sqrt(pow(2 * MT_PI, p_k) * p_det);
//   double p_normm = 1. / sqrt(pow(2 * MT_PI, p_k) * p_detm);

//   uint q_k = c.q_mean.N;
//   arr q_cov = c.q_cov;
//   // q_cov = eye(4);
//   arr q_covm = m * c.q_cov;
//   arr inv_q_cov = inverse(c.q_cov);
//   arr inv_q_covm = inverse(q_covm);
//   double q_det = determinant(c.q_cov);
//   double q_detm = determinant(q_covm);
//   double q_norm = 1. / sqrt(pow(2 * MT_PI, q_k) * q_det);
//   double q_normm = 1. / sqrt(pow(2 * MT_PI, q_k) * q_detm);

//   arr dpos(kflist.N, 3), dquat(kflist.N, 4);
//   arr dpos_centered(kflist.N, 3), dquat_centered(kflist.N, 4);
//   for(uint i = 0; i < kflist.N; i++) {
//     dpos[i]() = rec.query(STRING(s1 << "_dPos"), s2, kflist(i));
//     dquat[i]() = rec.query(STRING(s1 << "_dQuat"), s2, kflist(i));
//     dpos_centered[i]() = dpos[i]() - c.p_mean;
//     dquat_centered[i]() = dquat[i]() - c.q_mean;
//   }

//   double l = 0;
//   cout << "=== " << c.frame << " ========================================" << endl;
//   for(uint i = 0; i < kflist.N; i++) {
//     double p_mahal_dist_sq = scalarProduct(inv_p_cov, dpos_centered[i], dpos_centered[i]);
//     double q_mahal_dist_sq = scalarProduct(inv_q_cov, dquat_centered[i], dquat_centered[i]);
//     double p_mahal_dist = sqrt(p_mahal_dist_sq);
//     double q_mahal_dist = sqrt(q_mahal_dist_sq);


//     if(kflist(i) >= c.frame) {
//       cout << "AFTER" << endl;
//       // l += log(p_mahal_dist);
//       l += p_mahal_dist_sq;
//       // l += q_mahal_dist;
//       // l *= exp(-.5 * p_mahal_dist_sq) * p_norm;
//       // l *= exp(-.5 * q_mahal_dist_sq) * q_norm;
//       // l *= exp(-.5*scalarProduct(inv_q_cov, dquat_centered[i], dquat_centered[i])) * q_norm;
//     }
//     else {
//       cout << "BEFORE" << endl;
//       // l -= log(p_mahal_dist);
//       l -= p_mahal_dist_sq;
//       // l -= q_mahal_dist;
//       // Mahalanobis distance, i.e. value at 1std from the average?
//       // l *= sqrt(p_mahal_dist_sq);
//       // l *= sqrt(q_mahal_dist_sq);
//       // l *= exp(-.5*scalarProduct(inv_p_covm, dpos_centered[i], dpos_centered[i])) * p_normm;
//       // l *= exp(-.5*scalarProduct(inv_q_covm, dquat_centered[i], dquat_centered[i])) * q_normm;
//     }
//     cout << "p_mahal_dist: " << p_mahal_dist << endl;
//     cout << "p_mahal_dist_sq: " << p_mahal_dist_sq << endl;
//     cout << "--------" << endl;
//   }
//   cout << endl;
//   if(c.frame > 6000)
//     exit(0);
//   return l;
// }
// }}}

// KF_Model {{{
void KF_Model::setFeats(const MocapSeqL &seqlist) {
  watch::push();
  cout << "Setting features for seq no.";
  for(uint i = 0; i < seqlist.N; i++) {
    cout << " " << i << flush;
    setFeats(*seqlist(i));
  }
  cout << endl;
  cout << __PRETTY_FUNCTION__ << " done in " << watch::pop() << " sec." << endl;
}

MocapSeqL KF_Model::trainFilter(const MocapSeqL &seqlist) {
  MocapSeqL trainseqlist;
  for(MocapSeq *seq: seqlist) {
    seq->setAnn(target_str);
    if(seq->ann.N)
      trainseqlist.append(seq);
  }
  return trainseqlist;
}

void KF_Model::train(MocapData &mdata, const StringA &dsets) {
  MocapSeqL mseqlist;
  for(const String &dset: dsets) {
    MocapRec &mrec = mdata.rec(dset);
    mseqlist.append(mrec.seqlist());
    hash::append("params", std::string(dset));
  }

  // void KF_Model::train(const MocapSeqL &seqlist) {
  watch::push();
  if(!model().loadParams()) {

    MocapSeqL trainseqlist = trainFilter(mseqlist);
    cout << "num train sequences: " << trainseqlist.N << endl;

    setFeats(trainseqlist);

    // g4w.train(trainseqlist);
    // g4w.whiten(trainseqlist);

    watch::push();
    Model::DataL datal;
    for(MocapSeq *seq: trainseqlist) {
      Model::Data *data = model().newData();
      setData(*data, *seq);
      datal.append(data);
    }
    cout << "setting data done in " << watch::pop() << " sec." << endl;

    model().train(datal);

    listDelete(datal);
  }

  cout << __PRETTY_FUNCTION__ << " done in " << watch::pop() << " sec." << endl;
}

arr KF_Model::query(MocapSeq &seq, Thickness thickness) {
  setFeats(seq);
  // g4w.whiten(seq);

  Model::Data *data = model().newData();
  setData(*data, seq);

  arr y_thin;

  String &output = *params.get<String>("output");
  if(output == "maxlikelihood")
    y_thin = model().maxlikelihood(*data);
  else if(output == "marginal")
    y_thin = model().marginal(*data);
  else
    HALT("Output '" << output << "' not programmed.");
  delete data;

  if(thickness == THIN)
    return y_thin;
  return thicken(seq.nframes, y_thin);
}

double KF_Model::negloglikelihood(MocapSeq &seq) {
  setFeats(seq);
  // g4w.whiten(seq);

  Model::Data *data = model().newData();
  setData(*data, seq);
  double nll = model().negloglikelihood(*data);
  delete data;
  return nll;
}
// }}}
// KF_TS_Model {{{
KF_TS_Model::KF_TS_Model(): KF_Model(Target::TS) {}
KF_TS_Model::~KF_TS_Model() {}

void KF_TS_Model::setFeats(MocapSeq &seq) {
  StringA feats = MT::getParameter<StringA>("feats");

  seq.clearFeatData();
  for(const String &feat: feats)
    if(feat == "var")
      seq.appendVarPast();
    else if(feat == "lincoeff")
      seq.appendLinCoeffPast(false);
    // else if(feat == "ddist")
    //   seq.appendDDist();
    // else if(feat == "idist")
    //   seq.appendIDist();
  seq.data = seq.featdata;
}

void KF_TS_Model::play(MocapRec &rec, setf_t setf) {
  double col_on[3][3] = { { 1, 0, 0 },
                          { 0, 1, 0 },
                          { 0, 0, 1 } };
  double col_off[3] = { 1, 1, 1 };
  uint ncolors = sizeof(col_on) / sizeof(col_on[0][0]);

  const StringA &limbs = rec.id().agent_limbs();
  const StringA &parts = rec.id().object_parts();

  uint nlimbs = limbs.N;
  uint nparts = parts.N;
  uint nframes = rec.numFrames(THICK);

  CHECK(nlimbs <= ncolors, STRING("Not enough colors for " << nlimbs << " limbs; Add more."))

  MT::Array<StringA> digitsof(nlimbs);
  for(uint li = 0; li < nlimbs; li++)
    digitsof(li) = rec.id().sensorsof(limbs(li));

  arr interactions(nlimbs, nparts, nframes);
  for(uint li = 0; li < nlimbs; li++)
    for(uint pi = 0; pi < nparts; pi++)
      interactions[li][pi]() = maxlikelihood(rec, limbs(li), parts(pi), THICK);

  MocapGui gui(rec);
  gui.init_custom() = [&, col_on, nlimbs, limbs, digitsof] (ors::KinematicWorld &kw) {
    for(uint li = 0; li < nlimbs; li++)
      for(const String &digit: digitsof(li))
        for(ors::Shape *sh: kw.getBodyByName(digit)->shapes)
          memcpy(sh->color, col_on[li], 3 * sizeof(double));
  };
  gui.update_custom() = [&, col_on, col_off, nlimbs, nparts, parts, interactions] (ors::KinematicWorld &kw, uint f) {
    for(const String &part: parts)
      for(ors::Shape *sh: kw.getBodyByName(part)->shapes)
        memcpy(sh->color, col_off, 3 * sizeof(double));

    for(uint pi = 0; pi < nparts; pi++) {
      for(uint li = 0; li < nlimbs; li++) {
        if(interactions(li, pi, f)) {
          for(ors::Shape *sh: kw.getBodyByName(parts(pi))->shapes)
            memcpy(sh->color, col_on[li], 3 * sizeof(double));
          break;
        }
      }
    }
  };
  gui.play();
}

void KF_TS_Model::playTarget(MocapRec &rec) {
  NIY;
  // KF_Gui gui(rec);
  // double col_on[3][3] = { { 1, 0, 0 },
  //                         { 0, 1, 0 },
  //                         { 0, 0, 1 } };
  // double col_off[3] = { 1, 1, 1 };
  // uint ncolors = sizeof(col_on) / sizeof(col_on[0][0]);

  // Item *i = rec.agent_targets.getItem(Target_to_str(target));
  // const StringA &agents = (i != nullptr)? *i->getValue<StringA>(): StringA();
  // // const StringA &objects = rec.object_targets.getItem(target_str)?
  // //                         *rec.object_targets.getValue<StringA>(target_str):
  // //                         StringA();

  // uint nagents = agents.N;
  // CHECK(nagents <= ncolors, STRING("Not enough colors for " << nagents << " agents; Add more."))

  // StringA objectsensors, objectshapes;
  // objectsensors = rec.id().object_parts();
  // // for(const String &obj: objects)
  // //   objectsensors.append(rec.id().sensorsof(obj));
  // for(const String &objsensor: objectsensors)
  //   objectshapes.append(STRING("sh" << objsensor));
  // uint nobjectsensors = objectsensors.N;

  // MT::Array<StringA> agentshapes(nagents);
  // for(uint ai = 0; ai < nagents; ai++)
  //   for(const String &agentsensor: rec.id().sensorsof(agents(ai)))
  //     agentshapes(ai).append(STRING("sh" << agentsensor));

  // arr interactions;
  // for(const String &agent: agents)
  //   for(const String &objsensor: objectsensors)
  //     interactions.append(rec.ann(target, agent, objsensor));
  // interactions.reshape(nagents, nobjectsensors, interactions.N / (nagents * nobjectsensors));

  // gui.init();
  // while(true) {
  //   gui.step();
  //   if(gui.done())
  //     break;

  //   for(ors::Shape *sh: gui.kw().shapes)
  //     memcpy(sh->color, col_off, 3*sizeof(double));

  //   for(const String sensor: rec.id().sensors())
  //     gui.kw().getShapeByName(STRING("sh"<<sensor))->color[3] = rec.query("poseObs", sensor, gui.fnum()).scalar();

  //   for(uint ai = 0; ai < nagents; ai++)
  //     for(const String &agentshape: agentshapes(ai))
  //       memcpy(gui.kw().getShapeByName(agentshape)->color, col_on[ai], 3 * sizeof(double));

  //   for(uint osi = 0; osi < nobjectsensors; osi++) {
  //     for(uint ai = 0; ai < nagents; ai++) {
  //       if(interactions(ai, osi, gui.fnum())) {
  //         memcpy(gui.kw().getShapeByName(objectshapes(osi))->color, col_on[ai], 3 * sizeof(double));
  //         break;
  //       }
  //     }
  //   }
  //   gui.update();
  // }
  // gui.finally();
}

void KF_TS_Model::entitylists(MocapRec &mrec, StringA &elist1, StringA &elist2) {
  elist1 = mrec.id().agent_limbs();
  elist2 = mrec.id().object_parts();
}
// }}}
// KF_CE_Model {{{
KF_CE_Model::KF_CE_Model(): KF_Model(Target::CE) {}
KF_CE_Model::~KF_CE_Model() {}

void KF_CE_Model::setFeats(MocapSeq &seq) {
  StringA feats = MT::getParameter<StringA>("feats");

  seq.clearFeatData();
  for(const String &feat: feats)
    if(feat == "var")
      seq.appendVarFuture();
    else if(feat == "lincoeff")
      seq.appendLinCoeffFuture(false);
    // else if(feat == "ddist")
    //   seq.appendDDist();
    // else if(feat == "idist")
    //   seq.appendIDist();
  seq.data = seq.featdata;
}

void KF_CE_Model::play(MocapRec &rec, setf_t setf) {
  double col_on[3] = { 1, 1, 0 };
  double col_off[3] = { 1, 1, 1 };
  uint ncolors = sizeof(col_on) / sizeof(col_on[0]);

  if(rec.id().objects().N != 1)
    NIY;

  const StringA &parts = rec.id().object_parts();
  uint nparts = parts.N;
  // CHECK(nobjects <= ncolors, STRING("Not enough colors for " << nparts << " parts; Add more."))

  uint ninteractions = nparts * (nparts+1) / 2;
  uint nframes = rec.numFrames(THICK);
  arr interactions(nparts, nframes);
  interactions.setZero();
  for(uint pi1 = 0; pi1 < nparts; pi1++) {
    for(uint pi2 = 0; pi2 < pi1; pi2++) {
      MocapSeq *seq = rec.seq(parts(pi1), parts(pi2));
      arr q = query(*seq);
      interactions[pi1]() += q;
      interactions[pi2]() += q;
    }
  }

  MocapGui gui(rec);
  gui.update_custom() = [&, col_on, col_off, nparts, parts, interactions] (ors::KinematicWorld &kw, uint f) {
    for(uint pi = 0; pi < nparts; pi++) {
      const double *col = interactions(pi, f)? col_on: col_off;
      for(ors::Shape *sh: kw.getBodyByName(parts(pi))->shapes)
        memcpy(sh->color, col, 3 * sizeof(double));
    }
  };
  gui.play();
}

void KF_CE_Model::playTarget(MocapRec &rec) {
  NIY;
  // KF_Gui gui(rec);
  // double col_on[3][3] = { { 1, 0, 0 },
  //                         { 0, 1, 0 },
  //                         { 0, 0, 1 } };
  // double col_off[3] = { 1, 1, 1 };
  // uint ncolors = sizeof(col_on) / sizeof(col_on[0][0]);

  // const StringA &objects = rec.id().objects();

  // uint nobjects = objects.N;
  // CHECK(nobjects <= ncolors, STRING("Not enough colors for " << nobjects << " objects; Add more."))

  // StringA objectsensors, objectshapes;
  // for(const String &obj: objects)
  //   objectsensors.append(rec.id().sensorsof(obj));
  // for(const String &objsensor: objectsensors)
  //   objectshapes.append(STRING("sh" << objsensor));
  // uint nobjectsensors = objectsensors.N;

  // uint ninteractions = nobjectsensors * (nobjectsensors+1) / 2;
  // uint nframes = rec.numFrames();
  // arr interactions(ninteractions, nframes);
  // interactions.setZero();
  // for(uint os1 = 0; os1 < nobjectsensors; os1++)
  //   for(uint os2 = 0; os2 < os1; os2++) {
  //     arr &&inter = rec.ann(target, objectsensors(os1), objectsensors(os2));
  //     if(inter.N)
  //       interactions[os1*(os1+1)/2 + os2]() = inter;
  //   }
  // interactions.reshape(ninteractions, interactions.N / ninteractions);

  // gui.init();
  // while(true) {
  //   gui.step();
  //   if(gui.done())
  //     break;

  //   for(ors::Shape *sh: gui.kw().shapes)
  //     memcpy(sh->color, col_off, 3*sizeof(double));

  //   for(const String sensor: rec.id().sensors())
  //     gui.kw().getShapeByName(STRING("sh"<<sensor))->color[3] = rec.query("poseObs", sensor, gui.fnum()).scalar();

  //   for(uint os1 = 0; os1 < nobjectsensors; os1++) {
  //     gui.kw().getShapeByName(objectshapes(os1))->color[3] = rec.query("poseObs", objectsensors(os1), gui.fnum()).scalar();
  //     for(uint os2 = 0; os2 < os1; os2++) {
  //       if(interactions(os1*(os1+1)/2 + os2, gui.fnum())) {
  //         memcpy(gui.kw().getShapeByName(objectshapes(os1))->color, col_on[0], 3 * sizeof(double));
  //         memcpy(gui.kw().getShapeByName(objectshapes(os2))->color, col_on[0], 3 * sizeof(double));
  //       }
  //     }
  //   }
  //   gui.update();
  // }
  // gui.finally();
}

void KF_CE_Model::entitylists(MocapRec &mrec, StringA &elist1, StringA &elist2) {
  elist1 = elist2 = mrec.id().object_parts();
}

KF_ConstraintL KF_CE_Model::all_constraints(MocapRec &rec) {
  StringA objectsensors;
  for(const String &object: rec.id().objects())
    objectsensors.append(rec.id().sensorsof(object));
  uint nobjectsensors = objectsensors.N;

  CHECK(nobjectsensors > 1, "Constraint detection requires at least 2 objects");

  uint nframes = rec.numFrames(THICK);
  arr interactions(nobjectsensors, nobjectsensors, nframes);
  interactions.setZero();
  for(uint os1 = 0; os1 < nobjectsensors; os1++)
    for(uint os2 = 0; os2 < nobjectsensors; os2++) {
      if(os1 != os2) {
        MocapSeq *seq = rec.seq(objectsensors(os1), objectsensors(os2));
        interactions[os1][os2]() = query(*seq);
        delete seq;
      }
    }
  boolA constrained(nobjectsensors), other_constrained(nobjectsensors, nobjectsensors);
  constrained = false;
  other_constrained = false;

  uintA fs(nobjectsensors);
  fs = -1u;

  arr new_interactions(nobjectsensors, nframes);
  for(uint os1 = 0; os1 < nobjectsensors; os1++) {
    for(uint os2 = 0; os2 < nobjectsensors; os2++) {
      if(os1 != os2) {
        MocapSeq *seq = rec.seq(objectsensors(os1), objectsensors(os2));
        arr pair_inter = query(*seq);
        delete seq;

        int f = pair_inter.findValue(1);
        if(f != -1) {
          new_interactions[os1]().subRange(f, -1) = 1;
          new_interactions[os2]().subRange(f, -1) = 1;
          fs(os1) = (fs(os1) == -1u)? f: MT::MIN(fs(os1), f);
          fs(os2) = (fs(os2) == -1u)? f: MT::MIN(fs(os2), f);
          constrained(os1) = true;
          constrained(os2) = true;
          other_constrained(os1, os2) = true;
        }
      }
    }
  }

  cout << "CONSTRAINTS DETECTED:" << endl;
  for(uint os = 0; os < nobjectsensors; os++)
    cout << objectsensors(os) << ": " << (constrained(os)? "x": "_") << endl;

  cout << "CONSTRAINTS DETECTED:" << endl;
  for(uint os1 = 0; os1 < nobjectsensors; os1++) {
    cout << objectsensors(os1) << ":";
    for(uint os2 = 0; os2 < nobjectsensors; os2++)
      cout << (other_constrained(os1, os2)? " x": " _");
    cout << endl;
  }

  // uint root_os = rec.id().root_id();
  // if(root_os == -1u)
  //   for(uint os = 0; os < nobjectsensors; os++)
  //     if(fs(os) != -1u && (root_os == -1u || fs(os) < fs(root_os)))
  //       root_os = os;

  // CHECK(root_os != -1u, "No root detected?!");

  // cout << "root: " << objectsensors(root_os) << endl;

  // Actually computing constraints
  KF_ConstraintL clist;
  for(uint os1 = 0; os1 < nobjectsensors; os1++) {
    for(uint os2 = 0; os2 < nobjectsensors; os2++) {
      if(os1 != os2) {
        MocapSeq *seq = rec.seq(objectsensors(os1), objectsensors(os2));
        arr pair_inter = query(*seq);
        delete seq;
        int f = pair_inter.findValue(1);
        if(f != -1) {
          KF_Constraint *c = new KF_Constraint();

          c->frame = f;
          c->name = STRING("c" << os1 << "." << os2);
          c->main_obj = objectsensors(os1);
          c->obj = objectsensors(os2);

          rec.computeDPos(objectsensors(os1), objectsensors(os2));
          rec.computeDQuat(objectsensors(os1), objectsensors(os2));

          // TODO only take the observable stuff into account!!!!!
          arr dpos = rec.query(STRING(objectsensors(os1) << "_dPos"), objectsensors(os2)).subRange(f, -1);
          arr dquat = rec.query(STRING(objectsensors(os1) << "_dQuat"), objectsensors(os2)).subRange(f, -1);

          arr dposObs = rec.query(STRING(objectsensors(os1) << "_dPosObs"), objectsensors(os2)).subRange(f, -1);
          arr dquatObs = rec.query(STRING(objectsensors(os1) << "_dQuatObs"), objectsensors(os2)).subRange(f, -1);

          uint sum_dposObs = sum(dposObs);
          uint sum_dquatObs = sum(dquatObs);

          if(sum_dposObs == 0 || sum_dquatObs == 0) {
            cout << "ERROR: constraint between " << objectsensors(os1) << " and " << objectsensors(os2) << " detected, but can't extract transformation!" << endl;
            continue;
          }

          c->p_mean.resize(3); c->q_mean.resize(4);
          c->p_mean.setZero(); c->q_mean.setZero();
          for(uint t = 0; t < dpos.d0; t++) {
            if(dposObs(t))
              c->p_mean += dpos[t];
            if(dquatObs(t))
              c->q_mean += dquat[t];
          }
          c->p_mean /= sum(dposObs);
          c->q_mean /= sum(dquatObs);
          c->q_mean /= length(c->q_mean);

          // c->p_mean = sum(dpos, 0) / (double)dpos.d0;
          // // c->q_mean = sum(dquat, 0) / (double)dquat.d0;
          // // c->q_mean /= length(c->q_mean);
          // c->q_mean = sum(dquat, 0) / length(dquat);
          // c->q_mean /= length(c->q_mean);

          // cout << "constraint with " << objectsensors(os) << ":" << endl;
          // cout << " * dposObs.N: " << dposObs.N << endl;
          // cout << " * dquatObs.N: " << dquatObs.N << endl;
          // cout << " * nframes: " << rec.numFrames() << endl;
          // cout << " * fs: " << fs(os) << endl;
          // cout << " * sum1: " << sum(dposObs) << endl;
          // cout << " * sum2: " << sum(dquatObs) << endl;
          // cout << " * pos: " << c->p_mean << endl;
          // cout << " * quat: " << c->q_mean << endl;

          clist.append(c);
        }
      }
    }
  }
  return clist;
}

KF_ConstraintL KF_CE_Model::constraints(MocapRec &rec) {
  StringA objectsensors;
  for(const String &object: rec.id().objects())
    objectsensors.append(rec.id().sensorsof(object));
  uint nobjectsensors = objectsensors.N;

  CHECK(nobjectsensors > 1, "Constraint detection requires at least 2 objects");

  uint ninteractions = nobjectsensors * (nobjectsensors+1) / 2;
  uint nframes = rec.numFrames(THICK);
  arr interactions(ninteractions, nframes);
  interactions.setZero();
  for(uint os1 = 0; os1 < nobjectsensors; os1++)
    for(uint os2 = 0; os2 < os1; os2++) {
      MocapSeq *seq = rec.seq(objectsensors(os1), objectsensors(os2));
      interactions[os1*(os1+1)/2 + os2]() = query(*seq);
      delete seq;
    }
  interactions.reshape(ninteractions, interactions.N / ninteractions);

  boolA constrained(nobjectsensors), other_constrained(nobjectsensors, nobjectsensors);
  constrained = false;
  other_constrained = false;

  uintA fs(nobjectsensors);
  fs = -1u;

  arr new_interactions(nobjectsensors, nframes);
  for(uint os1 = 0; os1 < nobjectsensors; os1++) {
    for(uint os2 = 0; os2 < os1; os2++) {
      MocapSeq *seq = rec.seq(objectsensors(os1), objectsensors(os2));
      arr pair_inter = query(*seq);
      delete seq;

      int f = pair_inter.findValue(1);
      if(f != -1) {
        new_interactions[os1]().subRange(f, -1) = 1;
        new_interactions[os2]().subRange(f, -1) = 1;
        fs(os1) = (fs(os1) == -1u)? f: MT::MIN(fs(os1), f);
        fs(os2) = (fs(os2) == -1u)? f: MT::MIN(fs(os2), f);
        constrained(os1) = true;
        constrained(os2) = true;
        other_constrained(os1, os2) = true;
        other_constrained(os2, os1) = true;
      }
    }
  }

  cout << "CONSTRAINTS DETECTED:" << endl;
  for(uint os = 0; os < nobjectsensors; os++)
    cout << objectsensors(os) << ": " << (constrained(os)? "x": "_") << endl;

  cout << "CONSTRAINTS DETECTED:" << endl;
  for(uint os1 = 0; os1 < nobjectsensors; os1++) {
    cout << objectsensors(os1) << ":";
    for(uint os2 = 0; os2 < nobjectsensors; os2++)
      cout << (other_constrained(os1, os2)? " x": " _");
    cout << endl;
  }

  uint root_os = rec.id().root_id();
  if(root_os == -1u)
    for(uint os = 0; os < nobjectsensors; os++)
      if(fs(os) != -1u && (root_os == -1u || fs(os) < fs(root_os)))
        root_os = os;

  CHECK(root_os != -1u, "No root detected?!");

  cout << "root: " << objectsensors(root_os) << endl;

  // Actually computing constraints
  KF_ConstraintL clist;
  for(uint os = 0; os < nobjectsensors; os++) {
    if(os == root_os)
      continue;
    if(fs(os) == -1u) {
      cout << "WARNING: constraint with " << objectsensors(os) << " not detected!" << endl;
      continue;
    }

    KF_Constraint *c = new KF_Constraint();
    c->frame = fs(os);
    c->name = STRING("c" << os);
    c->main_obj = objectsensors(root_os);
    c->obj = objectsensors(os);

    rec.computeDPos(objectsensors(root_os), objectsensors(os));
    rec.computeDQuat(objectsensors(root_os), objectsensors(os));

    // TODO only take the observable stuff into account!!!!!
    arr dpos = rec.query(STRING(objectsensors(root_os) << "_dPos"), objectsensors(os)).subRange(fs(os), -1);
    arr dquat = rec.query(STRING(objectsensors(root_os) << "_dQuat"), objectsensors(os)).subRange(fs(os), -1);

    arr dposObs = rec.query(STRING(objectsensors(root_os) << "_dPosObs"), objectsensors(os)).subRange(fs(os), -1);
    arr dquatObs = rec.query(STRING(objectsensors(root_os) << "_dQuatObs"), objectsensors(os)).subRange(fs(os), -1);

    uint sum_dposObs = sum(dposObs);
    uint sum_dquatObs = sum(dquatObs);

    if(sum_dposObs == 0 || sum_dquatObs == 0) {
      cout << "ERROR: constraint with " << objectsensors(os) << " detected, but can't extract transformation!" << endl;
      continue;
    }

    c->p_mean.resize(3); c->q_mean.resize(4);
    c->p_mean.setZero(); c->q_mean.setZero();
    for(uint t = 0; t < dpos.d0; t++) {
      if(dposObs(t))
        c->p_mean += dpos[t];
      if(dquatObs(t))
        c->q_mean += dquat[t];
    }
    c->p_mean /= sum(dposObs);
    c->q_mean /= sum(dquatObs);
    c->q_mean /= length(c->q_mean);

    // c->p_mean = sum(dpos, 0) / (double)dpos.d0;
    // // c->q_mean = sum(dquat, 0) / (double)dquat.d0;
    // // c->q_mean /= length(c->q_mean);
    // c->q_mean = sum(dquat, 0) / length(dquat);
    // c->q_mean /= length(c->q_mean);

    // cout << "constraint with " << objectsensors(os) << ":" << endl;
    // cout << " * dposObs.N: " << dposObs.N << endl;
    // cout << " * dquatObs.N: " << dquatObs.N << endl;
    // cout << " * nframes: " << rec.numFrames() << endl;
    // cout << " * fs: " << fs(os) << endl;
    // cout << " * sum1: " << sum(dposObs) << endl;
    // cout << " * sum2: " << sum(dquatObs) << endl;
    // cout << " * pos: " << c->p_mean << endl;
    // cout << " * quat: " << c->q_mean << endl;

    clist.append(c);
  }
  return clist;
}

void KF_CE_Model::playConstraints(MocapRec &rec, const KF_ConstraintL &clist) {
  MocapGui gui(rec);

  gui.init_custom() = [&, clist] (ors::KinematicWorld &kw) {
    for(const KF_Constraint *c: clist) {
      ors::Body *b_main = gui.kw().getBodyByName(c->main_obj);
      ors::Body *b = gui.kw().getBodyByName(c->obj);

      ors::Transformation T;
      T.pos.set(c->p_mean.p);
      T.rot.set(c->q_mean.p);

      b->X = b_main->X * T;
    }
    kw.calc_fwdPropagateShapeFrames();
    kw.gl().update(NULL, true);
  };
  gui.show();
}

void KF_CE_Model::showConstraints(MocapRec &rec, const KF_ConstraintL &clist) {
  NIY;
  // KF_Gui gui(rec);

  // gui.init();
  // for(const KF_Constraint *c: clist) {
  //   ors::Body *b_main = gui.kw().getBodyByName(c->main_obj);
  //   ors::Body *b = gui.kw().getBodyByName(c->obj);

  //   ors::Transformation T;
  //   T.pos.set(c->p_mean.p);
  //   T.rot.set(c->q_mean.p);

  //   b->X = b_main->X * T;
  // }

  // gui.kw().calc_fwdPropagateShapeFrames();
  // gui.kw().gl().update(NULL, true);
  // while(!gui.done());
  // gui.finally();
}
// }}}

