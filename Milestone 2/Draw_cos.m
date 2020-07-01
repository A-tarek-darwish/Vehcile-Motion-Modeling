%% DESCRIPTION
%
% This is a script to initialise the coordinate system.
%
%% OUTPUT
%
% Coordinate system is drawn in the animation.
%
%
% V1.0 | 03-May-2016 | Martin Lankers | creation
%% 1.) Definitions
%% 1.) -General
x_ar = 0.5*(x_t_max_limit);              % Define x length of the arrow of the coordinate system
clr_cos = 'k';                           % Color of the mass
variable_cos = 'x';                      % Define the variable which is displayed at the coordinate system

%% 3.) Plot
%% 3.) -Draw cos
plotcos(x_ar,variable_cos,clr_cos,lnwdth,fntsz)
