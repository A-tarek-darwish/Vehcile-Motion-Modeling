%% DESCRIPTION
%
% This is a Script to solve the differential equation of a single mass 
% system.
%
%% OUTPUT
%
% Formatted figure of the displacement of a single mass system and its
% animation.
%
%
%
% V1.0 | 03-May-2016 | Martin Lankers | creation
%
%% Program
clear                                             % Delete Workspace
clc                                               % Clear Command Window
close all                                         % Close all figures

%% 1.) Definitions
%% 1.) -Parameter definition
mass                       = 750;                 % Mass of the body [kg]
stiffness                  = 50000;               % Stiffness Coefficient of spring [N/m]
damping                    = 0;                   % Damping coefficient of damper [Ns/m]
time                       = 0:0.01:1;            % Time [s]

x_0                        = 0.01;                % Initial Condition displacement
x_dot_0                    = 0;                   % Initial Condition velocity
    
%% 2.) Computing
w0 = [x_0 , x_dot_0];                             % Create a vector with initial conditions
A =  [0, 1;  (-1)*stiffness / mass, (-1)*damping / mass ]; % Create system Matrix
dw = @(t,w) A*w;                                  % Define derivative

%% 2.) -Numerical solution of the motion
[tsim,wsim] = ode45(dw,time,w0); % Calling numerical solver
            
