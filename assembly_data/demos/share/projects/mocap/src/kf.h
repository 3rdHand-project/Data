#pragma once
#include <Core/array.h>
#include <Core/keyValueGraph.h>
#include <Core/util.h>
// #include <Perception/g4data.h>
#include <Mocap/mocapdata.h>
#include <Mocap/utilAB.h>
#include <iomanip>
#include "optim.h"
#include "model.h"
// #include "kfgui.h"

#define VOTE_MODE_LIST(f) \
  f(MAJORITY) \
  f(UNANIMOUS) \
  f(NEGA_UNANIMOUS) \
  f(NO_VOTE_MODE)

// ENUM_MACRO(VoteMode, VOTEMODE_LIST);
#define VOTE_MODE_ENUM(elem) elem,
enum VoteMode { VOTE_MODE_LIST(VOTE_MODE_ENUM) };

const char *votemode_to_str(VoteMode mode);
VoteMode str_to_votemode(const char *str);

typedef std::function<uintA(const arr &)> setf_t;

auto onset = [](const arr &inter) {
  uintA kflist;
  for(uint i = 1; i < inter.d0; i++)
    if(inter(i) && !inter(i-1))
      kflist.append(i);
  return kflist;
};

auto offset = [](const arr &inter) {
  uintA kflist;
  for(uint i = 1; i < inter.d0; i++)
    if(!inter(i) && inter(i-1))
      kflist.append(i);
  return kflist;
};

auto bothset = [](const arr &inter) {
  uintA kflist;
  for(uint i = 1; i < inter.d0; i++)
    if(inter(i) != inter(i-1))
      kflist.append(i);
  return kflist;
};

auto noset = [](const arr &inter) {
  uintA kflist;
  return kflist;
};

setf_t polyset(setf_t monoset);

struct KF_Score;
struct KF_Constraint;
struct KF_Model;
struct KF_Label;

typedef MT::Array<KF_Constraint*>  KF_ConstraintL;
typedef MT::Array<KF_ConstraintL*> KF_ConstraintLL;

// KF_Score {{{
struct KF_Score {
  // old
  uint tot, totpos, totneg;
  uint scores[2][2];

  // new
  uint tp, fp, fn, tp_ss;

  KF_Score();
  KF_Score(const arr &t, const arr &y);

  void clear();

  double precision() const;
  double recall() const;
  double std() const;
  // double tpr() const;
  // double fpr() const;
  // double accuracy() const;
  // double f1() const;

  KF_Score operator+(const KF_Score &x);
  KF_Score& operator+=(const KF_Score &x);
  void write(std::ostream &os = std::cout) const;
};
stdOutPipe(KF_Score)
// KF_Score }}}
// ROCCurve {{{
// struct ROCCurve {
//     arr thresholds;
//     uintA ind;
//     MT::Array<KF_Score> scores;

//     void append(double thresh, const KF_Score &score);
//     double area();
//     void show(const char *pdf = NULL);
//   private:
//     void sort();
// };
// }}}
// KF_Constraint {{{
struct KF_Constraint {
  String name, main_obj, obj;

  uint frame;
  arr p_mean;
  arr q_mean;
  arr p_cov;
  arr q_cov;

  KF_Constraint() {}
  KF_Constraint(const KF_Constraint &c);
  ~KF_Constraint() {}

  void write(std::ostream &os = std::cout) const;
};
stdOutPipe(KF_Constraint)
// }}}
// KF_Model {{{
struct KF_Model: Parametric {
  public:
    Model *mod;
    const Target target;
    const String target_str;

    KF_Model(Target target);
    virtual ~KF_Model();

    Model &model();
    void setData(Model::Data &data, MocapSeq &seq);
    void setFeats(const MocapSeqL &seqlist);
    virtual void setFeats(MocapSeq &seq) = 0;

    void setDefaultParams();
    MocapSeqL trainFilter(const MocapSeqL &seqlist);
    // void train(const MocapSeqL &seqlist);
    void train(MocapData &mdata, const StringA &dsets);
    arr query(MocapSeq &seq, Thickness thickness = THICK);

    static arr thinnen(double nframes_thin, const arr &thick);
    static arr thicken(double nframes, const arr &thin);
    arr maxlikelihood(MocapRec &mrec, const String &limb, const String &part, Thickness thickness);
    arr marginal(MocapRec &mrec, const String &limb, const String &part);
    arr marginal2(MocapRec &mrec, const String &limb, const String &part);
    arr entropy(MocapRec &mrec, const String &limb, const String &part);

    virtual void entitylists(MocapRec &mrec, StringA &elist1, StringA &elist2) = 0;

    uint nseq_ground_truth(MocapRec &mrec, const String &limb, const String &part);
    uint nseq_maxlikelihood(MocapRec &mrec, const String &limb, const String &part);
    double nseq_expected(MocapRec &mrec, const String &limb, const String &part);
    uintA lseq_ground_truth(MocapRec &mrec, const String &limb, const String &part);
    uintA lseq_maxlikelihood(MocapRec &mrec, const String &limb, const String &part);
    uintA lseq_expected(MocapRec &mrec, const String &limb, const String &part);
    arr sample(MocapRec &mrec, const String &limb, const String &part);
    double precision(MocapRec &mrec, const String &limb, const String &part);
    double recall(MocapRec &mrec, const String &limb, const String &part);
    double estimate_entropy(MocapRec &mrec, const String &limb, const String &part);

    void plot_profile(MocapRec &mrec, const String &limb, const String &part);
    void plot_samples(MocapRec &mrec, const String &limb, const String &part);

    virtual void play(MocapRec &rec, setf_t setf = polyset(noset)) = 0;
    // virtual void playTarget(MocapRec &rec) = 0;

    // KF_Score statistics(const arr &seq1, const arr &seq2);
    KF_Score statistics(MocapRec &rec);
    KF_Score statistics(MocapRec &rec, const String &limb, const String &part);
    void print_statistics(MocapRec &rec);

    // KF_Score test(const MocapSeqL &seqlist);
    // KF_Score crossval(const MocapSeqL &seqlist);
    // ROCCurve ROC(const MocapSeqL &seqlist, const arr &thresholds);

    double negloglikelihood(MocapSeqL &seqlist);
    virtual double negloglikelihood(MocapSeq &seq);

    // uintA extract_constraints(const StringA &objects);
};
// KF_Model }}}
// KF_TS_Model {{{
struct KF_TS_Model: KF_Model {
  KF_TS_Model();
  ~KF_TS_Model();
  void setFeats(MocapSeq &seq);

  void play(MocapRec &rec, setf_t setf = polyset(noset));
  void playTarget(MocapRec &rec);
  void entitylists(MocapRec &mrec, StringA &elist1, StringA &elist2);
};
// }}}
// KF_CE_Model {{{
struct KF_CE_Model: KF_Model {
  KF_CE_Model();
  ~KF_CE_Model();
  void setFeats(MocapSeq &seq);

  void play(MocapRec &rec, setf_t setf = polyset(noset));
  void playTarget(MocapRec &rec);
  void entitylists(MocapRec &mrec, StringA &elist1, StringA &elist2);

  KF_ConstraintL constraints(MocapRec &rec);
  KF_ConstraintL all_constraints(MocapRec &rec);
  void showConstraints(MocapRec &rec, const KF_ConstraintL &clist);
  void playConstraints(MocapRec &rec, const KF_ConstraintL &clist);
};
// }}}
// KF_CE_CRF {{{
// struct KF_CE_CRF {
//   KF_ConstraintL extract(MocapRec &rec, KFModel &model, const StringA &objects);
//   KF_Constraint *computeConstraint(MocapRec &rec, const char *s1, const char *s2, uint f);
//   // KF_ConstraintL findConstraints(MocapRec &rec, KFModel &model, const StringA &agents, const StringA &objects);
//   // , const StringA &agents, const StringA &objects);
//   // KF_Constraint findConstraint(MocapRec &rec, const char *s1, const char *s2, const uintA &kflist);
//   // KF_Constraint computeConstraint(MocapRec &rec, const char *s1, const char *s2, const uintA &kflist);
//   // double computeConstraintScore(MocapRec &rec, const char *s1, const char *s2, const uintA &kflist, KF_Constraint &c);
// };
// }}}

