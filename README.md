# 3rdHand Data Repository

This repository contains the motion capture data collected by MLR and Inria.

## datasets

The `datasets/` folder is structured as follows:

 * `datasets/<scenario>/` contains all the demonstrations which share a common scenario (e.g. the assembly of a chair, or the cooperative assembly of a toolbox).
 * `datasets/<scenario>/scene.json` describes general info which is common to all the demonstrations of that same scenario.
 * `datasets/<scenario>/<demonstration>/` contains an actual demonstration. Please refer to the `README.md` file in each folder.

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
