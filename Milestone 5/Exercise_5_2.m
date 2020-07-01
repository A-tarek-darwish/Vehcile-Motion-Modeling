%% DESCRIPTION
%
% This is a Script to calculate the amplitude frequency response function
% of a two mass system under different kinds of excitations.
%
%% OUTPUT
%
% Formatted figure of the frequency response functions.
%
%% 1.) Parameters
% masses and inertias
mass                      = 1000;                   % Mass of the body [kg]
inertia                   = 1000;                   % Inertia of the body [kg*m^2]

% Stiffness and damping values
stiffness_f               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_f                 = 100;                   % Damping coefficient of damper [Ns/m]
stiffness_r               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_r                 = 100;                   % Damping coefficient of damper [Ns/m]

% Lengths center of gravity to front and rear end
length_f                  = 2.5;                    % Distance of the right spring-damper to the center of mass [m]
length_r                  = 2.5;                    % Distance of the left spring-damper to the center of mass [m]

omega = 0:0.01:50;                                  % Angular frequency vector [rad/s]

force = 2000;                                  % Constant force value [N]
length_force = 1.6;            % Distance point of force attack to center of gravity
    

%% 2.) Computing

% Set up system matrices
M = [mass 0 ; 0 inertia]; % Mass matrix

% Damping matrix
K = [damping_r+damping_f, length_f*damping_f-length_r*damping_r;...
    length_f*damping_f-length_r*damping_r, length_r^2*damping_r+length_f^2*damping_f];

% Stiffness matrix
C = [stiffness_r+stiffness_f, length_f*stiffness_f-length_r*stiffness_r;...
    length_f*stiffness_f-length_r*stiffness_r, length_r^2*stiffness_r+length_f^2*stiffness_f];

% Real excitation vector (splitted in sine and cosine vector)
h_c = [force; length_force*force];
h_s = [0;0];
    
%% 2.) Computing
%% 2.) -Calculation of FRF

% Initialisation of result vector
x = zeros(2,size(omega,2));

for k = 1:size(omega, 2) % for every omega
       
    % calculate complex excitation vector
    h_star = 1/2*(h_c - 1i*h_s);
    
    % calculate complex frequency response matrix
    inv_freq_matrix_complex = C - omega(k)^2*M + 1i*omega(k)*K;    % inverse matrix
    F_star =inv(inv_freq_matrix_complex);                          % matrix inversion
    
    % Calculation of complex solution
    x_star = F_star * h_star;
    
    % transfer to real solution
    x_c = real(2*x_star);
    x_s = -imag(2*x_star);
    
    % Calculation if amplitude
    x(:,k) = sqrt(x_c.^2 + x_s.^2);
    
end

% Splitting of Result vector into DOFs
x_hat = x(1,:);
phi_hat = x(2,:)/pi*180;

%% 3.) Plot
%% 3.) -Initialise figures
Exercise_5_2_initialize_figures                       % Runs the script to initialise the figures

hold on
set(graph_plot(1),'Parent',axes_graph(1), 'XData', omega, 'YData', x_hat);      % Set new Values to displacement graph
set(graph_plot(2),'Parent',axes_graph(2), 'XData', omega, 'YData', phi_hat);      % Set new Values to velocity graph
    
