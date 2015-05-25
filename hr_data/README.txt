
This folder contains two files:
	. plot_interaction_trajectories.m
	. tskc.mat
    
plot_interaction_trajectories.m is a self contained script that loads 
the file tskc.mat which contains the collected human-robot collaboration data
at IAS TU-Darmstadt.

The data contained in tskc.mat is plotted as 3D Cartesian trajectories of
the human wrist (tracked by motion capture) and the robot end-effector 
(by forward kinematics from the joint encoders), for each of the three 
types of collaborative tasks.

There are three types of human-robot collaboration tasks.
Each of the tasks are contained in a cell:

tskc{1}: plate handover
tskc{2}: screw handover
tskc{3}: screwdriver handover

Each field are described as:

tskc{k}.task: 
  The name of the task

tskc{k}.data: 
  A cell of D elements, where D is the number of demonstrations and each
  element is a M x N matrix of training data.
  The M rows:   
      correspond to the number of time steps of the trajectory.
  The N columns are given by 3+7+3 DoFs:
      3 XYZ DoFs of the tracked human wrist, followed by the 
      7 joint encoder measurements of the robot (starting from the shoulder link)
      3 XYZ DoFs of the robot end-effector (found by using forward
      kinematics)


