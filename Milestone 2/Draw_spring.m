%% DESCRIPTION
%
% This is a script to initialise the spring.
%
%% OUTPUT
%
% Spring is drawn in the animation.
%
%
% V1.0 | 03-May-2016 | Martin Lankers | creation
%% 1.) Definitions
%% 1.) -General
spring_number_windings = 8;                      % Number of spring windings
spring_radius = 0.1;                             % Radius of spring radius
phi_max = 2*pi*spring_number_windings;           % Calculate the maximum angle of spring rotations
phi_s = 0:pi/50:phi_max;                         % Define a vector in order to calculate y and z position of the spring vertices
y_offset_s = 0.2;                                % Spring y-offset              
y_pos_spring = spring_radius * sin(phi_s) + y_offset_s; % Calculate y position of spring vertices
z_pos_spring = spring_radius * cos(phi_s);       % Calculate z position of spring vertices
spring_foot = position_g(1) - dimension_g(1)/2;  % Position of the spring foot
spring_head = x_t(1) + dimension_m(1)/2;         % Initial position of spring head

%% 3.) Plot
%% 3.) -Draw spring
x_pos_spring = phi_s/phi_max * (spring_head - spring_foot) + spring_foot;               % Calculate initial x value for spring
plot3(axes_ani,x_pos_spring,y_pos_spring,z_pos_spring,'b','linewidth',lnwdth)           % Use plot3 function to plot spring
