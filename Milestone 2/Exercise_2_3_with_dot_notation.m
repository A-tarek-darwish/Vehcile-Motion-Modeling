%% DESCRIPTION
%
% This is a script to initialize the figures for an animation and a graph.
%
%% OUTPUT
%
% Formatted figure which can be used for an animation.
%
%
%
% Solve_equation_of_motion_analytically
%% 1.) Definitions
%% 1.) -General
clr = [236/255 237/255 237/255];                    % Background Color grey
unts = 'normalized';                                % Units for dimensions to normalized
lnwdth = 2;                                         % Linewidth 2
fntsz = 22;                                         % Fontsize 22

%% 1.) -Positions, titles and labels
pos_fig = [0.01 0.1 0.98 0.8];                      % Position and dimension of figure
title_graph = 'Displacement and velocity vs. time'; % Title of graph
title_ani = 'Tied one mass system';                 % Title of animation
xlabel_ani = 'Displacement x [m]';                  % Name of x-axis of Animation
xlabel_graph = 'Time t [s]';                        % Name of x-axis of Graph
ylabel_left = 'Displacement x [m]';                 % Name of first y-axis of Graph
ylabel_right = 'Velocity v [m/s]';                  % Name of second y-axis of Graph

%% 3.) Plot
%% 3.) -Initialise figures
fig = figure('color',clr,'units',unts,'position',pos_fig);        % Create a blank figure
subplot(1,2,2);                                                   % Divide figure into two subplots and select second plot
ax = gca;
ax.FontSize = fntsz;
ax.Color = clr;
xl = xlabel(xlabel_ani);
tit = title(title_graph);
xlim([time(1) time(end)]);
LineWidth = lnwdth;
box on

yyaxis left
ax.YColor ='k';
yl_left = ylabel(ylabel_left);
ax.YDir = 'reverse';
x_t_max_limit = (max(abs(x_t))+0.05*max(x_t));                  % Get Maximum of x_t and add 5 percent
ylim([-x_t_max_limit,x_t_max_limit]);

yyaxis right
ax.YColor ='r';
yl_right = ylabel(ylabel_right);
ax.YDir = 'reverse';
v_t_max_limit = (max(abs(v_t))+0.05*max(v_t));                  % Get Maximum of v_t and add 5 percent
ylim([-v_t_max_limit,v_t_max_limit]);

subplot(1,2,1)                                                      % Create and select first plot of figure
axes_ani = gca;                                                     % Save current axes
axes_ani.FontSize = fntsz;                                          % Set Fontsize of Animation
axes_ani.Color = clr;                                               % Set background color
axes_ani.XDir = 'reverse';                                          % Invert y-axis
xlim([-3*x_t_max_limit,4*x_t_max_limit])                        % Set the Limits of x-axis Animation (dependent of amplitude)
ylim([-2 2])                                                    % Set the Limit of y-axis

title(title_ani,'fontsize',fntsz)                               % Set title of Animation
xlabel(xlabel_ani,'fontsize',fntsz)                             % Set label x-axis of Animation


