%% DESCRIPTION
%
% This is a script to plot a cube in Matlab.
%
%% OUTPUT
%
% Formatted figure with a cube.
%
%
%% 1.) Definitions
%% 1.) -General
dimension = [1 1 1];  % Dimension of the cube
position  = [0 1 0];  % Position of the cube
clr_c = 'k';  % Color of the cube

%% 1.) -Needed Coordinates
x_c(1) = -dimension(1)/2;                   % x-coordinate of the cube is half the length
x_c(2) = dimension(1)/2;                    % x-coordinate of the cube is half the length
x_c(3) = dimension(1)/2;                    % x-coordinate of the cube is half the length
x_c(4) = -dimension(1)/2;                   % x-coordinate of the cube is half the length
x_c(5) = -dimension(1)/2;                   % x-coordinate of the cube is half the length
x_c(6) = dimension(1)/2;                    % x-coordinate of the cube is half the length
x_c(7) = dimension(1)/2;                    % x-coordinate of the cube is half the length
x_c(8) = -dimension(1)/2;                   % x-coordinate of the cube is half the length

y_c(1) = dimension(2)/2;                    % y-coordinate of the cube is half the length
y_c(2) = dimension(2)/2;                    % y-coordinate of the cube is half the length
y_c(3) = -dimension(2)/2;                   % y-coordinate of the cube is half the length
y_c(4) = -dimension(2)/2;                   % y-coordinate of the cube is half the length
y_c(5) = dimension(2)/2;                    % y-coordinate of the cube is half the length
y_c(6) = dimension(2)/2;                    % y-coordinate of the cube is half the length
y_c(7) = -dimension(2)/2;                   % y-coordinate of the cube is half the length
y_c(8) = -dimension(2)/2;                   % y-coordinate of the cube is half the length

z_c(1) = -dimension(3)/2;                   % z-coordinate of the cube is half the length
z_c(2) = -dimension(3)/2;                   % z-coordinate of the cube is half the length
z_c(3) = -dimension(3)/2;                   % z-coordinate of the cube is half the length
z_c(4) = -dimension(3)/2;                   % z-coordinate of the cube is half the length
z_c(5) = dimension(3)/2;                    % z-coordinate of the cube is half the length
z_c(6) = dimension(3)/2;                    % z-coordinate of the cube is half the length
z_c(7) = dimension(3)/2;                    % z-coordinate of the cube is half the length
z_c(8) = dimension(3)/2;                    % z-coordinate of the cube is half the length

vertices_c = [x_c(1) y_c(1) z_c(1);x_c(2) y_c(2) z_c(2);x_c(3) y_c(3) z_c(3);x_c(4) y_c(4) z_c(4);...
              x_c(5) y_c(5) z_c(5);x_c(6) y_c(6) z_c(6);x_c(7) y_c(7) z_c(7);x_c(8) y_c(8) z_c(8)];   % Define the eight corners of the cube

faces_c = [4 3 2 1;2 3 7 6;3 4 8 7;1 2 6 5;7 8 5 6;1 4 8 5];               % Define the faces of the cube by four corner for each face

%% 2.) -Translate vertices
vertices_translated = [vertices_c(:,1)+position(1) vertices_c(:,2)+position(2) vertices_c(:,3)+position(3)];

%% 3.) Plot
patch('Vertices', vertices_translated, 'Faces', faces_c, 'FaceColor', clr_c,'EdgeColor', clr_c);

xlabel('x')
ylabel('y')
zlabel('z')
