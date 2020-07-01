%% DESCRIPTION
%
% This is a Script to solve the differential equation of a force excited two degree of freedom
% system with an analytical approach.
%
%% OUTPUT
%
% Displacement of a force excited two degree of freedom system.
%
%% 1.) Definitions
%% 1.) -Parameter definition
% Masses and inertias
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
force                     = 2000;                                  % Constant force value [N]
length_force              = 1.6;            % Distance point of force attack to center of gravity
omega                     = 10;                % Angular frequency of the excitation [1/s]

% Time and initial conditions
time = 0:0.005:20;                                   % Time [s]
x_0 = 0.1;                                         % Initial Condition displacement [m]
x_dot_0 = 0;                                        % Initial Condition velocity [m/s]
phi_0 = 0.0;                                        % Initial Condition angle [rad]
phi_dot_0 = 0;                                      % Initial Condition angle velocity [rad/s]
    
%% 2.) Computing
%% 2.) Solving
% Set up system matrices
M = [mass 0 ; 0 inertia]; % Mass matrix

% Damping matrix
K = [damping_r+damping_f, length_f*damping_f-length_r*damping_r;...
    length_f*damping_f-length_r*damping_r, length_r^2*damping_r+length_f^2*damping_f];

% Stiffness matrix
C = [stiffness_r+stiffness_f, length_f*stiffness_f-length_r*stiffness_r;...
    length_f*stiffness_f-length_r*stiffness_r, length_r^2*stiffness_r+length_f^2*stiffness_f];

%% Calculate analytical homogenous solution
% solve eigenvalue problem for system
[eigenvector,lambda, cond] = polyeig(C,K,M);

% set up constant matrices
C_matr = [eigenvector(1,1) eigenvector(1,2) eigenvector(1,3) eigenvector(1,4);...
    eigenvector(2,1) eigenvector(2,2) eigenvector(2,3) eigenvector(2,4);...
    lambda(1)*eigenvector(1,1) lambda(2)*eigenvector(1,2) lambda(3)*eigenvector(1,3) lambda(4)*eigenvector(1,4);...
    lambda(1)*eigenvector(2,1) lambda(2)*eigenvector(2,2) lambda(3)*eigenvector(2,3) lambda(4)*eigenvector(2,4)];

% Real excitation vector (splitted in sine and cosine vector)
h_c = [force; length_force*force];
h_s = [0;0];

% calculate complex excitation vector
h_star = 1/2*(h_c - 1i*h_s);

% calculate complex frequency response matrix
inv_freq_matrix_complex = C - omega^2*M + 1i*omega*K;    % inverse matrix
F_star =inv(inv_freq_matrix_complex);                          % matrix inversion

% Calculation of complex solution
x_star = F_star * h_star;

% create vector of initial conditions
init_cond = [x_0-x_star(1)-conj(x_star(1)); phi_0-x_star(2)-conj(x_star(2));x_dot_0-1i*omega*x_star(1)+1i*omega*conj(x_star(1)); phi_dot_0-1i*omega*x_star(2)+1i*omega*conj(x_star(2))];

% solve constant values equation system
C_vec = C_matr\init_cond;

% Calculate homogeneous solution
x_t_h = C_vec(1)*exp(lambda(1)*time).*eigenvector(1,1) + C_vec(2)*exp(lambda(2)*time).*eigenvector(1,2) + C_vec(3)*exp(lambda(3)*time).*eigenvector(1,3) +...
    C_vec(4)*exp(lambda(4)*time).*eigenvector(1,4);
phi_t_h = C_vec(1)*exp(lambda(1)*time).*eigenvector(2,1) + C_vec(2)*exp(lambda(2)*time).*eigenvector(2,2) + C_vec(3)*exp(lambda(3)*time).*eigenvector(2,3) +...
    C_vec(4)*exp(lambda(4)*time).*eigenvector(2,4);

% Calculate particular solution
x_t_p = x_star(1)*exp(1i*omega*time)+conj(x_star(1))*exp(-1i*omega*time);
phi_t_p = x_star(2)*exp(1i*omega*time)+conj(x_star(2))*exp(-1i*omega*time);

% Calculate overall solution
x_t = real(x_t_h) + x_t_p;             % Car body displacement
phi_t = real(phi_t_h) + phi_t_p;         % Car body rotation

figure
plot(time,real(x_t_h),time,x_t_p,time,x_t)
legend([{'homogeneous'},{'particular'},{'overall'}],'Location','northeastoutside');

figure
plot(time,real(phi_t_h),time,phi_t_p,time,phi_t)
legend([{'homogeneous'},{'particular'},{'overall'}],'Location','northeastoutside');

% State the solution which you need if you want to analyse the motion only
% for a long period of time
x_t_long = x_t_p;             
phi_t_long = phi_t_p;        
