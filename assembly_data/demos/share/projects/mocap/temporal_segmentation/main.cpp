#include <Core/util.h>
#include <Core/util_t.h>
#include <Mocap/mocapdata.h>
#include <Mocap/jsondata.h>
#include <Mocap/g4data.h>
#include <Mocap/optidata.h>
#include <Mocap/utilAB.h>
#include "kf.h"
#include "optim.h"
#include "thread.h"

#include <unistd.h>

int main(int argc, char **argv) {
  MT::openConfigFile("config");

  StringA bases = MT::getParameter<StringA>("bases");
  StringA dsets_train = MT::getParameter<StringA>("dsets_train");
  String dset = MT::getParameter<String>("dset");

  MocapData mdata;
  // register data types
  // TODO automate this in MocapData::MocapData()
  // or even uniform data formats
  // mdata.reclist.append(new JsonRec());
  mdata.reclist.append(new G4Rec());
  mdata.reclist.append(new OptiRec());

  mdata.base().append(bases);

  KF_TS_Model ts_model;
  ts_model.train(mdata, dsets_train);
  MocapRec &mrec_test = mdata.rec(dset);
  ts_model.play(mrec_test);
}
