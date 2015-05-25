This program uses the temporal segmentation model to segment human manipulation
based on the correlated motion between hands and objects.  The results are
visualized in our own own ORS -- Open Robos Simulator.

Notice that it only makes sense to run this demo on the datasets containing
motion of the hands (most likly, toolbox-assembly datasets; to make sure,
visualize the datasets first using the appropriate demo).

# config

You can select which demonstration to run by changing the parameters `base` and
`dset`.  Set the `base` parameter to the absolute path of the data/ directory
in the 3rdHand_data repository, and the `dset` parameter to the name of the
demonstration which you want to play.

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

