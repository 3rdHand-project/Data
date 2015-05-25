function [tskc] = plot_interaction_trajectories( )
% This is a self contained script that loads the file tskc.mat
% which contains the collected human-robot collaboration data
% at IAS TU-Darmstadt.
%
% The data contained in tskc.mat is plotted as 3D Cartesian trajectories of
% the human wrist (tracked by motion capture) and the robot end-effector 
% (by forward kinematics from the joint encoders), for each of the three 
% types of collaborative tasks.
%
% There are three types of human-robot collaboration tasks.
% Each of the tasks are contained in a cell:
%
% tskc{1}: plate handover
% tskc{2}: screw handover
% tskc{3}: screwdriver handover
% 
% Each field are described as:
% 
% tskc{k}.task: 
%   The name of the task
%
% tskc{k}.data: 
%   A cell of D elements, where D is the number of demonstrations and each
%   element is a M x N matrix of training data.
%   The M rows:   
%       correspond to the number of time steps of the trajectory.
%   The N columns are given by 3+7+3 DoFs:
%       3 XYZ DoFs of the tracked human wrist, followed by the 
%       7 joint encoder measurements of the robot (starting from the shoulder link)
%       3 XYZ DoFs of the robot end-effector (found by using forward
%       kinematics)
%
% The function 'myplot' receives the M x N matrix of training data and
% plots in 3D the first 3 columns representing the XYZ coordinates of the
% human with the last 3 columns representing the XYZ coordinates of the
% robot

    clear; close all; clc; dbstop if error;
    
    % load the raw trajectories of different human-robot interactions.
    load('tskc.mat');
    
    zheight = -0.41;
    sground = struct('LineWidth', 4, 'Color', [0.8 0.8 0.8], 'LineStyle', '-');
    
    % define linestyles
    sh1 = struct('LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
    sr1 = struct('LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
    sh2 = struct('LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
    sr2 = struct('LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
    sh3 = struct('LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
    sr3 = struct('LineWidth', 2, 'Color', 'k', 'LineStyle', '-');

    % 1st interaction pattern
    h.threeA = figure; hold on; grid on; axis 'equal'; box off; 
    title 'Plate handover';
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

    ir = 430:numel(tskc{1}.data{1}(:,1)); % show only the part with variance
    ir = 1:numel(tskc{1}.data{1}(:,1)); % show whole traj
    
    myplot(tskc{1}.data, sh1, sr2, sground, zheight, ir);
    view([0.5 -0.5 0.3]);        
    xlim([ 0.4   1.2]);
    ylim([-0.8 0.6]);
    zlim([zheight 0.4]);
    text(0.8, 0.2, 0.4, 'Human');
    text(0.8, -0.5, 0.4, 'Robot');    

    % 2nd interaction pattern
    h.threeB = figure; hold on; grid on; axis 'equal'; box off;
    title 'Screw handover';
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    view([0.5 0.5 0.1]);
    
    ir = 120:numel(tskc{1}.data{1}(:,1)); % show only the part with variance
    ir = 1:numel(tskc{1}.data{1}(:,1)); % show whole traj
    myplot(tskc{2}.data, sh1, sr2, sground, zheight, ir);
    view([0.5 -0.5 0.3]);        
    xlim([ 0.4   1.2]);
    ylim([-0.7 0.6]);
    zlim([zheight 0.2]);
    text(0.8, 0.2, 0.05, 'Human');
    text(0.8, -0.5, 0.2, 'Robot');     
    
    % 3rd interaction pattern
    h.threeC = figure; hold on; grid on; axis 'equal'; box off; 
    title 'Screwdriver handover';
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    view([0.5 0.5 0.1]);    
    
    ir = 1:numel(tskc{1}.data{1}(:,1));
    myplot(tskc{3}.data, sh1, sr2, sground, zheight, ir);
    view([0.5 -0.5 0.3]);        
    xlim([ 0.4   1.2]);
    ylim([-0.6 0.6]);
    zlim([zheight 0.0]);
    text(0.8, 0.2, -0.01, 'Human');
    text(0.8, -0.5, -0.01, 'Robot'); 
end

function [] = myplot(data, sh1, sr1, sground, zheight, ir)
% This function receives the M x N matrix of training data and
% plots in 3D the first 3 columns representing the XYZ coordinates of the
% human with the last 3 columns representing the XYZ coordinates of the
% robot

    % This is an auxiliar function used to plot the shadow of the trajectories.
    z = @(data) zheight * ones(numel(data),1);

    ndemo = numel(data);
    for d = 1:ndemo
        h_human_motion = plot3( data{d}(:,1),   data{d}(:,2),   data{d}(:,3));
        h_human_shadow = plot3( data{d}(:,1),   data{d}(:,2),   z(data{d}(:,2)));
        h_robot_motion = plot3( data{d}(ir,11), data{d}(ir,12), data{d}(ir,13));
        h_robot_shadow = plot3( data{d}(ir,11), data{d}(ir,12), z(data{d}(ir,13)));

        set(h_human_motion, sh1);
        set(h_human_shadow, sground);
        set(h_robot_motion, sr1);
        set(h_robot_shadow, sground);
    end
end
