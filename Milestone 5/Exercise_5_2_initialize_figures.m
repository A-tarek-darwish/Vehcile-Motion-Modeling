%% DESCRIPTION
%
% This is a script to initialize the figures for an animation and a graph.
%
%% OUTPUT
%
% Formatted figure which can be used for an animation.
%
%% 1.) Definitions
%% 1.) -General
clr = [236/255 237/255 237/255];                    % Background Color grey
unts = 'normalized';                                % Units for dimensions to normalized
lnwdth = 2;                                         % Linewidth 2
fntsz = 22;                                         % Fontsize 22

%% 1.) -Positions, titles and labels
pos_fig = [0.01 0.2 0.75 0.65];                      % Position and dimension of figure
title_graph = 'Frequency response'; % Title of graph
xlabel_graph = '\Omega [rad/s]';                            % x-axis of Graph: angular frequency of excitation
ylabel_graph{1} = 'Amplitude x_{hat} [m]';                      % Name of first y-axis of Graph
ylabel_graph{2} = 'Amplitude \phi_{hat} [°]';                 % Name of second y-axis of Graph

%% 3.) Plot
%% 3.) -Initialise figures
fig = figure('color',clr,'units',unts,'position',pos_fig);      % Create a blank figure
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

x_t_max_limit = (max(abs(x_hat))+0.05*max(x_hat));                  % Get Maximum of x_t and add 5 percent
ylim(axes_graph(1),[0,x_t_max_limit]);             % Limit first y-axis
% 
v_t_max_limit = (max(abs(phi_hat))+0.05*max(phi_hat));                  % Get Maximum of v_t and add 5 percent
ylim(axes_graph(2),[0,v_t_max_limit]);             % Limit second y-axis
% 
xlim(axes_graph(1),[omega(1) omega(end)]);                        % Limit the time axis
