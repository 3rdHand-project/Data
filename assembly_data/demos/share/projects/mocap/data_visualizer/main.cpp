#include <Core/util.h>
#include <Core/util_t.h>
#include <Mocap/mocapdata.h>
#include <Mocap/jsondata.h>
#include <Mocap/g4data.h>
#include <Mocap/optidata.h>
#include <Mocap/mocapgui.h>

int main(int argc, char **argv) {
  MT::openConfigFile("config");

  StringA bases = MT::getParameter<StringA>("bases");
  String dset = MT::getParameter<String>("dset");
  String target = MT::getParameter<String>("target");

  cout << "bases:   " << bases << endl;
  cout << "dataset: " << dset << endl;
  cout << "target:  " << target << endl;

  MocapData mdata;
  
  // register data types
  // TODO automate this in MocapData::MocapData()
  // or even uniform data formats
  // mdata.reclist.append(new JsonRec());
  mdata.reclist.append(new G4Rec());
  mdata.reclist.append(new OptiRec());

  // set up data base directory & load dataset
  mdata.base().append(bases);
  MocapRec &mrec = mdata.rec(dset);

  // play data
  MocapGui mgui(mrec);
  mgui.play(str_to_Target(target));
}
