%% DESCRIPTION
%
% This is a script to draw the mass.
%
%% OUTPUT
%
% Mass is drawn in the animation.
%
%
% V1.0 | 03-May-2016 | Martin Lankers | creation
%% 1.) Definitions
%% 1.) -General
dimension_m = [1.2*x_t_max_limit 1 0.2];                        % Length, width and height of the mass
position_m = [x_t(1) 0 0];                                      % Initial position of the mass
clr_m = 'r';                                                    % Color of the mass

%% 3.) Plot
%% 3.) -Draw mass
plotcube(axes_ani,dimension_m,position_m,clr_m)                 % Initialise the mass
