%% DESCRIPTION
%
% This is a Script to solve the differential equation of a two degree of freedom
% system using the ODE funtion.
%
%% OUTPUT
%
% Displacement of a a two degree of freedom system.
%
%% 1.) Definitions
%% 1.) -Parameter definition
% Masses and inertias
mass                      = 1000;                 % Mass of the body [kg]
inertia                   = 1000;                 % Inertia of the body [kg*m^2]

% Stiffness and damping values
stiffness_f               = 60000;                % Stiffness coefficient of spring [N/m]
damping_f                 = 0;                    % Damping coefficient of damper [Ns/m]
stiffness_r               = 60000;                % Stiffness coefficient of spring [N/m]
damping_r                 = 0;                    % Damping coefficient of damper [Ns/m]

% Lengths center of gravity to front and rear end
length_f                  = 2.5;                  % Distance of the right spring-damper to the center of mass [m]
length_r                  = 2.5;                  % Distance of the left spring-damper to the center of mass [m]

% Time and initial conditions
time = 0:0.005:1;                                 % Time [s]
x_0 = 0.1;                                        % Initial Condition displacement [m]
x_dot_0 = 2;                                      % Initial Condition velocity [m/s]
phi_0 = 0.1;                                      % Initial Condition angle [rad]
phi_dot_0 = 1;                                    % Initial Condition angle velocity [rad/s]
    
%% 2.) Computing
%% 2.) Solving
%% 2.) Computing
%% 2.) -Numerical solution of the motion
w0 = [x_0,phi_0,x_dot_0,phi_dot_0];                               % Create a vector with initial conditions

A = [ 0 0 1 0; 0 0 0 1;...
    ((-stiffness_r - stiffness_f)/mass) ((-stiffness_f*length_f + stiffness_r*length_r)/mass) ((-damping_r - damping_f)/mass) ((-damping_f*length_f + damping_r*length_r)/mass);...
    ((-stiffness_f*length_f + stiffness_r*length_r)/inertia) ((-stiffness_f*length_f^2 - stiffness_r*length_r^2)/inertia) ((-damping_f*length_f + damping_r*length_r)/inertia) ((-damping_f*length_f^2 - damping_r*length_r^2)/inertia)];

% Return state space representation
dw = @(t,w) A*w;

[tsim,wsim] = ode45(dw,time,w0); % Calling numerical solver
