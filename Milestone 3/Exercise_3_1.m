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
damping                    = 1000;                % Damping coefficient of damper [Ns/m]
time                       = 0:0.01:15;           % Time [s]

x_0                        = 0.0;                 % Initial Condition displacement
x_dot_0                    = 0;                   % Initial Condition velocity

mass_extruded              = 0.03;                % Mass of the excitation [kg]
omega                      = 6.28;                % Angular frequency of the excitation [1/s]
radius                     = 0.24;                % Radius of the excitation [m]

force = omega^2*radius*mass_extruded;             % Calculate the force with given parameters
    
%% 2.) Computing
%% 2.) -Parameter calculation
dampingcoefficient = damping/(2*mass);             % Calculate damping coefficient
angulareigenfrequency = sqrt(stiffness/mass);      % Calculate angular eigenfrequency

%% 2.) -Calculation of the characteristic polynomial
lambda = roots([1, 2*dampingcoefficient, angulareigenfrequency^2]); %Calculate roots

%% 2.) -Calculation of the constants
x_max = sqrt(force^2/((damping*omega)^2 + (stiffness-mass*omega^2)^2));   % maximum amplitude of particular solution
phi = atan((damping*omega)/(stiffness-mass*omega^2));                     % phase shift of the particular solution

k1 = x_0 - (lambda(1)*x_0 - x_dot_0 - x_max*omega*sin(-phi) - x_max*lambda(1)*cos(-phi)) / (lambda(1) - lambda(2))- x_max*cos(-phi);      % Calculation of the constant 1
k2 = (lambda(1)*x_0 - x_dot_0 - x_max*omega*sin(-phi) - x_max*lambda(1)*cos(-phi)) / (lambda(1) - lambda(2));                             % Calculation of the constant 2

%% 2.) -Calculation of the solution
x_t_h = k1 * exp(lambda(1)*time) + k2 * exp(lambda(2)*time);                         % % Calculate homogeneous solution
v_t_h = k1 * lambda(1) * exp(lambda(1)*time) + k2 * lambda(2) * exp(lambda(2)*time); % % Calculate homogeneous solution

x_t_p = x_max*cos(omega*time-phi);                                          % Calculate particular solution
v_t_p = -x_max*omega*sin(omega*time-phi);                                   % Calculate particular solution

x_t = real(x_t_h) + x_t_p;            % Calculate overall solution by superposition of homogeneous and particular solution
v_t = real(v_t_h) + v_t_p;            % Calculate overall solution by superposition of homogeneous and particular solution

plot(time,x_t_h,time,x_t_p,time,x_t)
