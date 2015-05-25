#include <Core/util.h>
#include <Core/util_t.h>
#include <Mocap/mocapdata.h>
#include <Mocap/jsondata.h>
#include <Mocap/optidata.h>
#include <Mocap/g4data.h>
#include <Mocap/mocapgui.h>
#include "kf.h"
// #include "gui.h"
#include "optim.h"
// #include "util.h"
// #include "mocapplayer.h"
// #include "kfgui.h"
#include "thread.h"

#include <unistd.h>

void save_json(const KF_ConstraintLL &cllist);

int main(int argc, char **argv) {
  String::readStopSymbols = " \n\r";
  MT::initCmdLine(argc, argv);
  String config = MT::getParameter<String>("config", STRING("config"));
  MT::openConfigFile(config);

  StringA bases = MT::getParameter<StringA>("bases");
  StringA dsets_train = MT::getParameter<StringA>("dsets_train");
  StringA dsets = MT::getParameter<StringA>("dsets");

  MocapData mdata;
  // register data types
  // TODO automate this in MocapData::MocapData()
  // or even uniform data formats
  // mdata.reclist.append(new JsonRec());
  mdata.reclist.append(new G4Rec());
  mdata.reclist.append(new OptiRec());

  mdata.base().append(bases);

  // MocapRecL mreclist_train;
  // MocapSeqL mseqlist_train;
  // for(const String &dset_train: dsets_train) {
  //   MocapRec &mrec_train = mdata.rec(dset_train);
  //   mreclist_train.append(&mrec_train);
  //   mseqlist_train.append(mrec_train.seqlist());
  // }

  // cout << "num train recordings: " << mreclist_train.N << endl;
  // cout << "num train sequences:  " << mseqlist_train.N << endl;

  KF_CE_Model ce_model;
  // ce_model.train(mseqlist_train);
  ce_model.train(mdata, dsets_train);

  String task;
  MT::getParameter<String>("task", task);

  // KF_ConstraintLL cllist;
  for(const String &dset: dsets) {
    MocapRec &mrec = mdata.rec(dset);
    const KF_ConstraintL &clist = ce_model.constraints(mrec);
    // cout << "CONSTRAINTS" << endl;
    // if(task == "play")
    //   ce_model.play(mrec);
    // else if(task == "show") {
    // cout << "SHOW" << endl;
      ce_model.playConstraints(mrec, clist);
    // }
    // else {
    //   cout << "NONE" << endl;
    // }
    // else if(task == "show")
    //   ce_model.showConstraints(mrec, clist);
    // cllist.append(new KF_ConstraintL(clist));
  }
  // if(task == "json")
  //   save_json(cllist);
}

// void save_json(const KF_ConstraintLL &cllist) {
//   Json::Value json_root, json_demo, json_constraint, json_array;

//   json_root.clear();
//   for(const KF_ConstraintL *clist: cllist) {
//     json_demo.clear();
//     for(const KF_Constraint *c: *clist) {
//       cout << "== == ==" << endl;
//       cout << c->main_obj << endl;
//       cout << c->obj << endl;
//       cout << c->p_mean << endl;
//       cout << c->q_mean << endl;
//       cout << "== == ==" << endl;

//       json_constraint.clear();
//       json_constraint["object_reference"] = (const char *)c->main_obj;
//       json_constraint["object_relative"] = (const char *)c->obj;
//       json_constraint["time"] = c->frame;

//       json_array.clear();
//       json_array.append(c->p_mean(0));
//       json_array.append(c->p_mean(1));
//       json_array.append(c->p_mean(2));
//       json_constraint["constraint"].append(json_array);

//       json_array.clear();
//       json_array.append(c->q_mean(1));
//       json_array.append(c->q_mean(2));
//       json_array.append(c->q_mean(3));
//       json_array.append(c->q_mean(0));
//       json_constraint["constraint"].append(json_array);;

//       json_demo.append(json_constraint);
//     }
//     json_root.append(json_demo);
//   }

//   String json_filename = MT::getParameter<String>("json_filename");

//   ofstream fout;
//   MT::open(fout, json_filename);

//   fout << json_root;
//   // Json::FastWriter writer;
//   // fout << writer.write(root);

//   fout.close();
// }
