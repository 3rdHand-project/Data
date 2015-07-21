# Robot Programming from Demonstration, Feedback and Transfer (IROS 2015)

This folder contains source data used to compute the results presented in the paper *Robot Programming from Demonstration, Feedback and Transfer*, IROS 2015.

## Experimental protocol
We recruited 14 subjects and evaluate what is the impact of each process into the overall efficiency of the system from the high-level knowledge acquisition to its execution in a MoveIt simulated environment of the Baxter robot.
Four different exercises, labeled from A to D, are considered and they all consist into programming a robot to assemble either a chair or a bench. In all exercises the subject has to understand the environment around the robot and select a plan that avoids obstacles and workspace limits of the robot. The four exercises are listed below:

* **Chair cold start:** In this exercise the subject builds the chair by starting directly with the *Task Refining* process from no demonstration. He needs to create all constraints and plan from scratch using only the GUI. Our hypothesis is that this process is not intuitive, will take a longer time and require more effort from the user.
* **Chair bootstrap:** This exercise uses the *Task Learning* process to compute constraints and assembly plan from demonstrations and then bootstrap the *Task Refining* process. The subject is asked to perform corrections on the data learned from his own demonstrations in order to assemble the chair successfully in simulation. Our hypothesis is that the initial demonstration process is very intuitive and will give a good quality overall plan, but will not have enough precision nor consider the runtime environment to fulfill the task. So the second process of refining will give this extra advantage.
* **Bench cold start:** Similarly to the first exercise, this one consists into assembling the bench from no demonstration, from scratch through the GUI only.
* **Bench bootstrap:** This exercise illustrates knowledge transfer from the chair to the bench. The subject is asked to assemble the bench using the plan of the chair that he has already corrected. We preload the plan corresponding to the most precise constraints between exercises A or B. Thus the only corrections that he needs to provide correspond to the difference between the chair and the bench that cannot be computed automatically. Our hypothesis is that using the demonstrations of a different object allows nevertheless to save time and effort compared to a cold start.

The subject starts by recording 3 demonstrations of the chair assembly acquired with an Optitrack tracking system. The second step of the protocol consists into training the subject to the GUI during 20 minutes with synthetic example constraints and plans. These data contain all possible cases to train into all features of the GUI: create and delete a constraint, reorder the plan, correct precision, and also get familiar with 3D navigation ensured with a 6D space mouse.

Then the subject performs randomly either exercises A and then B or B and then A. Both exercises must lead to a successful assembly of the chair in simulation. In case of failure, the user performs new corrections and retries as many times as needed. Finally the subject performs randomly either exercises C and then D or D and then C, both exercises leading to a successful assembly of the bench.

## Content description

The file **all_subjects_dump.yaml** is a dump of the data recorded by the system while the subject was performing the exercises, i.e. the time and effort in terms of number of clicks while correcting the constraints though the GUI and simulating the assembly. All the 14 subjects are aggregated into this YAML file. YAML format is easily readable in Python using the module of the same name.

## Content format

TODO
