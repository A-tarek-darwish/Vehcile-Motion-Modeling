%% DESCRIPTION
%
% This is a script to draw the ground.
%
%% OUTPUT
%
% Ground is drawn in the animation.
%
%
% V1.0 | 03-May-2016 | Martin Lankers | creation
%% 1.) Definitions
%% 1.) -General
dimension_g = [0.25*x_t_max_limit 2 2];                         % Length, width and height of the ground
position_g = [abs(min(x_t))*2.5 0 0];                           % Position of the ground depending on the minimum displacement of the mass
clr_g = 'k';                                                    % Color of the ground

%% 3.) Plot
%% 3.) -Draw ground
plotcube(axes_ani,dimension_g,position_g,clr_g)                 % Initialise the ground
