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
%% Program
clear                                             % Delete Workspace
clc                                               % Clear Command Window
close all                                         % Close all figures

%% 1.) Definitions
%% 1.) -Parameter definition
mass                       = 750;                 % Mass of the body [kg]
stiffness                  = 50000;               % Stiffness Coefficient of spring [N/m]
damping                    = 0;                   % Damping coefficient of damper [Ns/m]
time                       = 0:0.01:2;            % Time [s]

x_0                        = 0.01;                % Initial Condition displacement
x_dot_0                    = 0.2;                   % Initial Condition velocity
    
%% 1.) -Symbolic function definition
t = sym('t');                                     % Define time variable
x = symfun(sym(str2sym('x(t)')), [t]);            % Define dependent variable
Dx = diff(x,1);                                   % Define first derivation
D2x = diff(x,2);                                  % Define second derivation

%% 2.) Computing
%% 2.) -Solve the equation
x = dsolve( mass*D2x + damping*Dx + stiffness*x == 0  ,x(0)==x_0  ,Dx(0) == x_dot_0,  't');  % Initial conditions units are [m] and [m/s] respectively
        
%% 2.) -Evaluate the equation
x_fun = matlabFunction(x);                        % Create function handle for x
x_dot_fun = matlabFunction(diff(x));              % Create function handle for x_dot

x_t = x_fun(time);                                % Evaluate function at points "time"
v_t = x_dot_fun(time);                            % Evaluate function at points "time"
              
%% 2.) -Calculate amplitude
x_roof = max(abs(x_t));                           % Find the maximum of the absolute values of the displacement

%% 2.) -Calculate time period
[maxima,max_location] = findpeaks(x_t,time);      % Use findpeaks to find local maxima
time_period = diff(max_location);                 % Subtract the locations of maxima

%% 2.) -Calculate frequency
frequency = 1/time_period(1);                     % Frequency is the reciprocal of the time period

%% 2.) -Calculate angular eigenfrequency
angular_eigenfrequency = 2*pi*frequency;          % Use constant pi to calculate angular eigenfrequency

%% 2.) -Calculate phase angle
temp_variable = diff(sign(v_t));                  % Find differences in the sign
indx_up = find(temp_variable>0);                  % Find zero crossings from negative to positive
indx_down = find(temp_variable<0);                % Find zero crossings from positive to negative
indx_up_time = time(indx_up);                     % Find corresponding time 
indx_down_time = time(indx_down);                 % Find corresponding time 
first_zero_crossing = min([indx_up_time indx_down_time]);       % Calculate the first zero crossing
phase_angle_rad = angular_eigenfrequency*first_zero_crossing;   % Calculate the phase angle in radian
phase_angle_degree = 180/pi*phase_angle_rad;                    % Convert phase angle in degree
