%% DESCRIPTION
%
% This is a script to initialize the figures for an animation and a graph.
%
%% OUTPUT
%
% Formatted figure which can be used for an animation.
%
%
Solve_equation_of_motion_analytically
%% 1.) Definitions
%% 1.) -General
clr = [236/255 237/255 237/255];                    % Background Color grey
unts = 'normalized';                                % Units for dimensions to normalized
lnwdth = 2;                                         % Linewidth 2
fntsz = 22;                                         % Fontsize 22

%% 1.) -Positions, titles and labels
pos_fig = [0.01 0.1 0.98 0.8];                      % Position and dimension of figure
title_graph = 'Displacement and velocity vs. time'; % Title of graph
title_ani = 'Tied one mass system';        % Title of animation
xlabel_ani = 'Displacement x [m]';                  % Name of x-axis of Animation
xlabel_graph = 'Time t [s]';                        % Name of x-axis of Graph
ylabel_graph{1} = 'Displacement x [m]';             % Name of first y-axis of Graph
ylabel_graph{2} = 'Velocity v [m/s]';               % Name of second y-axis of Graph

%% 3.) Plot
%% 3.) -Initialise figures
fig = figure('color',clr,'units',unts,'position',pos_fig);      % Create a blank figure
subplot(1,2,2)                                                  % Divide figure into two subplots and select second plot
graph_plot = plot(1,1,1,1);                                     % Initialise graph at second position
set(graph_plot(1),'color', 'k','linewidth',lnwdth);             % Set Color and linewidth of first plot
set(gca, 'Color', clr)
set(graph_plot(2),'color', 'r','linewidth',lnwdth);             % Set Color and linewidth of second plot
axes_graph(1) = gca;                                            % Save first yaxis
set(axes_graph(1),'FontSize',fntsz);                            % Set Fontsize of x and y-axes
axes_graph(2) = axes('Position',axes_graph(1).Position,...
    'YAxisLocation','right','YColor','r','Color','none','XTickLabel',[],'fontsize',fntsz);   % Create and save second yaxis
linkaxes([axes_graph(1) axes_graph(2) ],'x');                   % Link both x Axes to each other
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz)             % Label x-axis
ylabel(axes_graph(1),ylabel_graph{1},'fontsize',fntsz)          % Label y-axis
ylabel(axes_graph(2),ylabel_graph{2},'fontsize',fntsz)          % Label y-axis
title(title_graph,'fontsize',fntsz);                            % Title of Graph
set(axes_graph(1),'Ydir','reverse')                             % Invert y-axis
set(axes_graph(2),'Ydir','reverse')                             % Invert y-axis

x_t_max_limit = (max(abs(x_t))+0.05*max(x_t));                  % Get Maximum of x_t and add 5 percent
ylim(axes_graph(1),[-x_t_max_limit,x_t_max_limit]);             % Limit first y-axis

v_t_max_limit = (max(abs(v_t))+0.05*max(v_t));                  % Get Maximum of v_t and add 5 percent
ylim(axes_graph(2),[-v_t_max_limit,v_t_max_limit]);             % Limit second y-axis

xlim(axes_graph(1),[time(1) time(end)]);                        % Limit the time axis

subplot(1,2,1)                                                  % Create and select first plot of figure
axes_ani = gca;                                                 % Save current axes
set(axes_ani,'FontSize',fntsz);                                 % Set Fontsize of Animation
set(gca, 'Color', clr)                                          % Set background color
set(axes_ani,'Xdir','reverse')                                  % Invert y-axis
xlim([-3*x_t_max_limit,4*x_t_max_limit])                        % Set the Limits of x-axis Animation (dependent of amplitude)
ylim([-2 2])                                                    % Set the Limit of y-axis

title(title_ani,'fontsize',fntsz)                               % Set title of Animation
xlabel(xlabel_ani,'fontsize',fntsz)                             % Set label x-axis of Animation
