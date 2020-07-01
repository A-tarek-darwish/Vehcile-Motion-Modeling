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
%% Program
clear                                             % Delete Workspace
clc                                               % Clear Command Window
close all                                         % Close all figures

%% 1.) Definitions
%% 1.) -Parameter definition
mass                       = 750;                 % Mass of the body [kg]
stiffness                  = 50000;               % Stiffness Coefficient of spring [N/m]
damping1                   = 0;                   % Damping coefficient of damper [Ns/m]
damping2                   = 1000;                   % Damping coefficient of damper [Ns/m]
time                       = 0:0.01:1;            % Time [s]

x_0                        = 0.01;                % Initial Condition displacement
x_dot_0                    = 0;                   % Initial Condition velocity

y_0                        = 0.01;                % Initial Condition displacement
y_dot_0                    = 0;                   % Initial Condition velocity
    
%% 1.) -Symbolic function definition
t = sym('t');                                     % Define time variable
x = symfun(sym(str2sym('x(t)')), [t]);            % Define dependent variable
Dx = diff(x,1);                                   % Define first derivation
D2x = diff(x,2);                                  % Define second derivation

t = sym('t');                                     % Define time variable
y = symfun(sym(str2sym('y(t)')), [t]);            % Define dependent variable
Dy = diff(y,1);                                   % Define first derivation
D2y = diff(y,2);                                  % Define second derivation

%% 2.) Computing
%% 2.) -Solve the equation
x_res = dsolve( mass*D2x + damping1*Dx + stiffness*x == 0  ,x(0)==x_0  ,Dx(0) == x_dot_0,  't');  % Initial conditions units are [m] and [m/s] respectively
y_res = dsolve( mass*D2y + damping2*Dy + stiffness*y == 0  ,y(0)==y_0  ,Dy(0) == y_dot_0,  't');  % Initial conditions units are [m] and [m/s] respectively
        
%% 2.) -Evaluate the equation
x_fun = matlabFunction(x_res);                        % Create function handle for x
x_dot_fun = matlabFunction(diff(x_res));              % Create function handle for x_dot

x_t = x_fun(time);                                % Evaluate function at points "time"
vx_t = x_dot_fun(time);                            % Evaluate function at points "time"

y_fun = matlabFunction(y_res);                        % Create function handle for x
y_dot_fun = matlabFunction(diff(y_res));              % Create function handle for x_dot

y_t = y_fun(time);                                % Evaluate function at points "time"
vy_t = y_dot_fun(time);                            % Evaluate function at points "time"

difference = x_t-y_t;
