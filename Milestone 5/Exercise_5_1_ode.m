%% DESCRIPTION
%
% This is a Script to solve the differential equation of a force excited two degree of freedom
% system using the ODE funtion.
%
%% OUTPUT
%
% Displacement of a force excited two degree of freedom system.
%
%% 1.) Definitions
%% 1.) -Parameter definition
% Masses and inertias
mass                      = 1000;                   % Mass of the body [kg]
inertia                   = 1000;                   % Inertia of the body [kg*m^2]

% Stiffness and damping values
stiffness_f               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_f                 = 100;                   % Damping coefficient of damper [Ns/m]
stiffness_r               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_r                 = 100;                   % Damping coefficient of damper [Ns/m]

% Lengths center of gravity to front and rear end
length_f                  = 2.5;                    % Distance of the right spring-damper to the center of mass [m]
length_r                  = 2.5;                    % Distance of the left spring-damper to the center of mass [m]
force = 2000;                                  % Constant force value [N]
length_force = 1.6;            % Distance point of force attack to center of gravity
omega                      = 10;                % Angular frequency of the excitation [1/s]

% Time and initial conditions
time = 0:0.005:20;                                   % Time [s]
x_0 = 0.1;                                         % Initial Condition displacement [m]
x_dot_0 = 0;                                        % Initial Condition velocity [m/s]
phi_0 = 0.0;                                        % Initial Condition angle [rad]
phi_dot_0 = 0;                                      % Initial Condition angle velocity [rad/s]
    
%% 2.) Computing
%% 2.) -Numerical solution of the motion
w0 = [x_0,phi_0,x_dot_0,phi_dot_0];                               % Create a vector with initial conditions
% Set up system matrix
A = [ 0 0 1 0; 0 0 0 1;...
    ((-stiffness_r - stiffness_f)/mass) ((-stiffness_f*length_f + stiffness_r*length_r)/mass) ((-damping_r - damping_f)/mass) ((-damping_f*length_f + damping_r*length_r)/mass);...
    ((-stiffness_f*length_f + stiffness_r*length_r)/inertia) ((-stiffness_f*length_f^2 - stiffness_r*length_r^2)/inertia) ((-damping_f*length_f + damping_r*length_r)/inertia) ((-damping_f*length_f^2 - damping_r*length_r^2)/inertia)];
% Set up excitation vector
B = [0; 0; (force/mass); (force*length_force/inertia)];

dw = @(t,w) A*w+B*cos(omega*t);

[tsim,wsim] = ode45(dw,time,w0);

figure
plot(time,wsim(:,1))
legend('overall','Location','northeastoutside');

figure
plot(time,wsim(:,2))
legend('overall','Location','northeastoutside');
