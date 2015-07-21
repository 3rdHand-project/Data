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

The file [all_subjects_dump.yaml](all_subjects_dump.yaml) is a dump of the data recorded by the system while the subject was performing the exercises, i.e. the time and effort in terms of number of clicks while correcting the constraints though the GUI and simulating the assembly. All the 14 subjects are aggregated into this YAML file. YAML format is easily readable in Python using the module of the same name.

## Content format

The file format is described through a moke example:

```
john:    # Unique subject name
  blue_chair:    # Object (blue_chair or green_bench, redundant with the exercise)
    'A':    # Exercise (A, B, C, or D)
      num_attempts: 2    # john solved exercise A in 2 attempts
      '0':  # Assembly attempt #0 over 2
        effort:    # Measurement of effort provided to correct the constraints (number of clicks and time)
          actions: {kfviews: 26,    # The user visualized some assembly steps (keyframes) 26 times 
                    reordering: 4}    # The user made 4 changes in the plan by reordering its steps
          back: {click: 1,    # The back has been clicked once (to give it the focus or trigger the contextual menu)
                 kfviews: 6,    # Housekeeping data
                 movex: 5,    # The back has been moved 5 times along X
                 movey: 8,    # The back has been moved 8 times along Y
                 movez: 6,    # The back has been moved 6 times along Z
                 parent: 1,    # The parent of the back has been changed once
                 reordering: 4,    # Housekeeping data
                 rotx: 1,    # The back has been rotated once around X
                 roty: 1,    # The back has been rotated once around Y
                 rotz: 1}    # The back has been rotated once around Z
          seating: {click: 6}    # The seating has been clicked 6 times (and no operation was made on its constraints)
        corrected_plan:    # The assembly plan for this assembly attempt #0 after correction
        - [mc2, seating, leg1]    # Step 1 attributes constraint mc2 to leg1 with parent seating 
        - [mc3, seating, leg2]    # Step 2 attributes constraint mc3 to leg2 with parent seating 
        - [mc4, seating, leg3]    # Step 3 attributes constraint mc4 to leg3 with parent seating 
        - [mc5, seating, leg4]    # Step 4 attributes constraint mc5 to leg4 with parent seating 
        - [mc1, seating, back]    # Step 5 attributes constraint mc1 to back with parent seating 
        - wait    # The last step of all plans is always wait
        corrected_precision:    # The precision measured for this assembly attempt #0 after correction with respect to ground truth
          orientation: [0.04632476885626439, 0.04892203013021544, 0.03541035763741666,
            0.009341378860458691, 0.002442223469800918]
          position: [0.007897460392127359, 0.010246950765959597, 0.0067082042990198225,
            0.005385164216364273, 0.003352116981704301]
          reference_id: [c3, c1, c2, c4, c0]    # Housekeeping data
        corrected_st:    # Symbols table giving geometrical values to each constraint ID
          mc1: [-0.0014894702471792698,-0.003, 0.11986521631479263, -0.00019683582650031894,
            -0.00604905653744936, 0.0011892152251675725, 0.9999809861183167, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]    # position (x, y, z in metres), orientation (the quaternion x, y, z, w), degrees of freedom in translation (-x, +x, -y, +y, -z, +z in metres) and degrees of freedom in orientation (RPY -x, +x, -y, +y, -z, +z in radians)
          mc2:    # Next constraint in this symbols table...
        loadings: [1424278244.2105181]    # Timestamps of all clicks on the "Load constraint" button loading either an empty plan or a plan extracted from demonstrations depending on the exercise
        writings: [1424279232.2969055, 1424279306.7999418]  # Timestamps of all clicks on the "Submit constraints" button giving the signal that corrections are finished for this attempt. When length(writings)>1, the last timestamp corresponds to a validation of DoFs from the experimentor and must be ignored.
      plan:    # Housekeeping data
      symbols:    # Housekeeping data
      state: success    # After simulation attempt #1 was a success 
      '1':  # Assembly attempt #1 over 2 (i.e. the last one)
        # [...]
    'B':    # Next exercise...
```
