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
%close all                                         % Close all figures

%% 1.) Definitions
%% 1.) -Parameter definition
mass                       = 750;                 % Mass of the body [kg]
stiffness                  = 50000;               % Stiffness Coefficient of spring [N/m]
damping                    = 0;                   % Damping coefficient of damper [Ns/m]
time                       = 0:0.01:1;            % Time [s]

x_0                        = 0.01;                % Initial Condition displacement
x_dot_0                    = 0;                   % Initial Condition velocity
    
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
              
%% 3.) Plot
%% 3.) -Plot parameters
clr = [236/255 237/255 237/255];                    % Background Color grey
unts = 'normalized';                                % Units for dimensions to normalized
lnwdth = 2;                                         % Linewidth 2
fntsz = 22;                                         % Fontsize 18
pos_fig = [0.01 0.1 0.98 0.8];                      % Position and dimension of figure
title_graph = 'Displacement and velocity vs. time'; % Title of graph
xlabel_graph = 'Time t [s]';                        % Name of x-axis of Graph
ylabel_graph{1} = 'Displacement x [m]';             % Name of first y-axis of Graph
ylabel_graph{2} = 'Velocity v [m/s]';               % Name of second y-axis of Graph

%% 3.) -Plot graph
fig = figure('color',clr,'units',unts,'position',pos_fig);              % Create a blank figure
[axes_graph, Line1,Line2] = plotyy(time,x_t,time,v_t);                  % Plot figure with two y axes
set(Line1,'color', 'k','linewidth',lnwdth);                             % Set Color and linewidth of first plot
set(Line2,'color', 'r','linewidth',lnwdth);                             % Set Color and linewidth of second plot
set(axes_graph(1),'Ycolor', 'k','linewidth',lnwdth,'fontsize',fntsz);   % Set Color and linewidth of first plot
set(axes_graph(2),'Ycolor', 'r','linewidth',lnwdth,'fontsize',fntsz);   % Set Color and linewidth of second plot
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);                    % Label x-axis
ylabel(axes_graph(1),ylabel_graph{1},'fontsize',fntsz);                 % Label y-axis
ylabel(axes_graph(2),ylabel_graph{2},'fontsize',fntsz);                 % Label y-axis
title(title_graph);                                                     % Title of the graph 

x_t_max_limit = (max(abs(x_t))+0.05*max(abs(x_t)));                     % Get Maximum of x_t and add 5 percent
ylim(axes_graph(1),[-x_t_max_limit,x_t_max_limit]);                     % Limit first y-axis
v_t_max_limit = (max(abs(v_t))+0.05*max(abs(v_t)));                     % Get Maximum of v_t and add 5 percent
ylim(axes_graph(2),[-v_t_max_limit,v_t_max_limit]);                     % Limit second y-axis
