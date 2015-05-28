# 3rdHand Data Repository

This repository contains the motion capture data collected by the 3rdHand consortium.

## datasets

The `datasets/` folder is structured as follows:

 * `datasets/<scenario>/` contains all the demonstrations which share a common scenario (e.g. the assembly of a chair, or the cooperative assembly of a toolbox).
 * `datasets/<scenario>/scene.json` describes general info which is common to all the demonstrations of that same scenario.
 * `datasets/<scenario>/<demonstration>/` contains an actual demonstration. Please refer to the `README.md` file in each folder.

### scene.json

The scene.json file contains general info about a set of demonstrations.
Generically speaking, this includes the objects and agents involved in a
demonstration.

The json data structure contains two fields: `objects` and `agents`.
The `objects` field contains an entry for each object, which can be indexed by the object's name

import json

f = open("scene.json")
data = json.load(f)

data["objects"] # contains an entry for each object-related sensor.

Each entry is indexed by an object's name which is composed by the
concatenation of the hierarchical components of the object, e.g. "/chair/leg",
"/toolbox/side_long", "/box", etc..

Each entry contains the following fields:

data["object"]["<obj_name>"]["type"] # This defines the type of an object, and
allows to identify whether different objects are in fact interchangeable.

data["object"]["<obj_name>"]["mesh"] # Is the name of the mesh-file for the
object.

data["object"]["<obj_name>"]["transformation"] # Defines the transformation
between the mesh's coordinate system and the sensor's position.

data["agents"] # contains an entry for each agent-related sensor.

As above, each entry is indexed by the agent sensor's identity which, in the
case of agent sensors, is divided into <agent_name>, <limb_name> and
<finger_name>, e.g. "/human/rh/index", "/human/lh/thumb", etc..

Because there are no mesh models for the fingers of the human, each entry only
contains the following field:

data["agents"]["/<agent>/<limb>/<digit>"]["type"] # This defines the type of
the sensor. In the case of a human-related sensor, this is always "digit". The
field exists exclusively to provide a list of the existing human-related
sensors, and potentially for future development when more attributes will be
associated to the sensors.

## meshes

The `meshes/` folder contains object models, which are used to play the
demonstrations back in our software.

The `.ors` format is one used internally by the MLR group to describe kinematic
chains, joints, and degrees of freedom. You most likely won't be interested in
those models, since they are very coarse and correspond to minor entities (e.g.
the finger tips).

More likely, you might be interested in using the `.stl` files, which are
proper mesh files for the various parts of the objects of interest (e.g. the
toolbox, the chair, and the stool).

## demos

The `demos/` folder contains software to visualize the dataset, and demos which
showcase the research that has been developed for the 3rdHand project.

Please refer to the `README.md` file you will find there.
