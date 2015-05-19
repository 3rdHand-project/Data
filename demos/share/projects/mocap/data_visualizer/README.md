This program plays the motion capture data collected by MLR and Inria using our
own ORS -- Open Robos Simulator.

# config

You can select which demonstration to run by changing the parameters `base` and
`dset`.  Set the `base` parameter to the absolute path of the data/ directory
in the 3rdHand_data repository, and the `dset` parameter to the name of the
demonstration which you want to play.

The `target` variable allows you to visualize the provided labels for a given
task, where `TS` stands for the temporal_segmentation task, and `CE` stands for
the constraint_extraction task.  In the case where a demonstration does not
contain the labels corresponding to a specified task (e.g. demonstrations where
the hands are not tracked can not be used for the temporal_segmentation task),
the program will exit with a message indicating this event.

In the config file, if you set the same variable more than once, only the first
will be taken into account, and you can't comment assignments out.  If you
have multiple assignments you want to use at different times, just move the one
which you want to use at any given time on top of the others.

# gui

The gui which appears when you run ./x.exe can be easily controlled using mouse
and keyboard.

 * mouse left click and drag changes the viewing angle.
 * mouse right click changes the focus of the window.
 * mouse wheel controls the zoom.
 * key 0 inverts the direction of play (forward or backwards).
 * keys 1-9 determine the speed of reproduction.
 * keys j,k manually move backward and forward through the frames.
 * spacebar pauses and plays.
 * key q closes the program.

