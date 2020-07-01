%% DESCRIPTION
%
% This is a script to initialise the damper.
%
%% OUTPUT
%
% Damper is drawn in the animation.
%
%
% V1.0 | 03-May-2016 | Martin Lankers | creation
%% 1.) Definitions
%% 1.) -General
clr_d = 'k';                                                    % Color of the mass
y_offset_d = -0.2;                                              % Damper y-offset
damper_foot = position_g(1) - dimension_g(1)/2;                 % Position of the damper foot
damper_head = x_t(1) + dimension_m(1)/2;                        % Initial position of damper head
stroke_length_max = abs(min(x_t + dimension_m(1)/2)) + damper_foot;           % Maximum stroke length

%% 3.) Plot
%% 3.) -Draw damper
plotdamper(stroke_length_max,damper_foot,damper_head,y_offset_d,clr_d,lnwdth) % Plot damper
