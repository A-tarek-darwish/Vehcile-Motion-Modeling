%% DESCRIPTION
%
% This is a Script to solve the differential equation of a force excited single mass 
% system.
%
%% OUTPUT
%
% Formatted figure of the displacement of a single mass system and its
% animation.
%
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

mass_extruded              = 0.03;                % Mass of the excitation [kg]
omega                      = 6.28;                % Angular frequency of the excitation [1/s]
radius                     = 0.24;                % Radius of the excitation [m]

force = omega^2*radius*mass_extruded;             % Calculate the force with given parameters
    
%% 2.) Computing
w0 = [x_0 , x_dot_0];                             % Create a vector with initial conditions
A = [0, 1;  (-1)*stiffness / mass, (-1)*damping / mass ]; % Create system Matrix
B = [0 ; force / mass];                           % Create excitation vector
dw = @(t,w) A*w + B*cos(omega*t);                 % Define derivative    

%% 2.) -Numerical solution of the motion
[tsim,wsim] = ode45(dw,time,w0); % Calling numerical solver        
