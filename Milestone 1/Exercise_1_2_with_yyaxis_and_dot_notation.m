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
time                       = 0:0.01:1;            % Time [s]

x_0                        = 0.01;                % Initial Condition displacement
x_dot_0                    = 0;                   % Initial Condition velocity
    
%% 1.) -Symbolic function definition
syms x(t)                                         % Define dependent variable
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
              
%% 3.) Plot
%% 3.) -Creating the Plot
clr = [236/255 237/255 237/255];                    % Background Color grey
unts = 'normalized';                                % Units for dimensions to normalized
lnwdth = 2;                                         % Linewidth 2
fntsz = 22;                                         % Fontsize 18
pos_fig = [0.01 0.1 0.98 0.8];                      % Position and dimension of figure

fig = figure('color',clr,'units',unts,'position',pos_fig);              % Create a blank figure
ax = gca;
ax.LineWidth = lnwdth;
ax.FontSize = fntsz;

yyaxis left                                      % Plot figure with two y axes
ax = gca;
ax.YColor = 'k';
x_plot = plot(time, x_t);

yyaxis right
ax = gca;
ax.YColor = 'r';
v_plot = plot(time, v_t);

xl = xlabel("Time t in [s]");                    % Name of x-axis of Graph
yyaxis left
yl_left = ylabel("Displacement x in [m]");       % Name of left y-axis of Graph
yyaxis right
yl_right = ylabel("Velocity v in [m/s]");        % Name of right y-axis of Graph

ti = title("Displacement and Velocity vs. Time");% Title of graph

%% 3.) -Plot formatting
LineWidth_pref = 2;
x_plot.LineWidth = LineWidth_pref;
v_plot.LineWidth = LineWidth_pref;
x_plot.Color = 'k';
v_plot.Color = 'r';

x_t_max_limit = max(abs(x_t)) + 0.05*max(abs(x_t));
v_t_max_limit = max(abs(v_t)) + 0.05*max(abs(v_t));

yyaxis left
ylim([-x_t_max_limit, x_t_max_limit])
yyaxis right
ylim([-v_t_max_limit, v_t_max_limit])
