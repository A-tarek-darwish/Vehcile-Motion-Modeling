%% DESCRIPTION
%
% This is a Script to solve the differential equation of a force excited two degree of freedom
% system with different approaches.
%
%% OUTPUT
%
% Displacement of a force excited two degree of freedom system.
%
%% 1.) Definitions
%% 1.) -Parameter definition
length_force_0 = 1.6;            % Initial distance point of force to center of gravity
length_force_lower_bound       = -4;  % Lower bound of force length
length_force_upper_bound       = 4; % Upper bound of force length

force                      = 2000;                % Force on the body [N]
omega                      = 2.4;                 % Angular frequency of the excitation [1/s]

length_force_vector = length_force_lower_bound:0.01:length_force_upper_bound; % Define a force vector from the lower to the upper bound with an inkrement of 0.01
cost = nan(size(length_force_vector));  % Initialize the cost vector
k = 0;  % Initialize the counting variable
for length_force  = length_force_lower_bound:0.01:length_force_upper_bound  % Loop over the length force
    k = k + 1;  % Counting variable goes one up
    cost(k) = Exercise_5_3_force_length_obj(length_force,omega);  % Excecute the objective function
end
plot(length_force_vector,cost)  % Plot and inspect the cost function. Plot the cost values over the length_force_vector
   
%% 2.) Computing
%% 2.) Solving with local search
f = @(x) Exercise_5_3_force_length_obj(x,omega);        % Define function handle for the local search
length_force_opt_one_omega = fminsearch(f,length_force_0);

%% 2.) Solving with local search for different omega
kk = 0;  % Initialise counting variable
for omega = 1:0.1:30   % Define a omega vector from 1 to 30 with an inkrement of 0.1 and loop over that vector
    f = @(x) Exercise_5_3_force_length_obj(x,omega);      % Define function handle for the local search
    kk = kk + 1;  % Counting variable goes one up
    % You could use fmincon on your local MATLAB version. It is not
    % supported on the edx gui. You can also define boundary conditions
    % with the fmincon command.
    %length_force_opt_fmincon(omega) = fmincon(f,length_force_0,[],[],[],[],length_force_lower_bound,length_force_upper_bound); 
    length_force_opt_fminsearch(k) = fminsearch(f,length_force_0);  % Use fminsearch to find the minimum of the cost function
end


function cost = Exercise_5_3_force_length_obj(length_force_0,omega)
% % DESCRIPTION
% 
% This is a Script to solve the differential equation of a two degree of freedom
% system with different approaches.
% 
% % OUTPUT
% 
% Displacement of a a two degree of freedom system.
% 
% % VERSION
%             author: Martin Lankers (Martin.Lankers.de)
%      creation date: 25-April-2018
%     Matlab version: 2017a
% 
% % REVISION
% 
% V1.0 | 25-April-2018 | Martin Lankers | creation

% Masses and inertias
mass                      = 1000;                   % Mass of the body [kg]
inertia                   = 1000;                   % Inertia of the body [kg*m^2]

% Stiffness and damping values
stiffness_f               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_f                 = 1000;                   % Damping coefficient of damper [Ns/m]
stiffness_r               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_r                 = 1000;                   % Damping coefficient of damper [Ns/m]

% Lengths center of gravity to front and rear end
length_f                  = 2.5;                    % Distance of the right spring-damper to the center of mass [m]
length_r                  = 1.5;                    % Distance of the left spring-damper to the center of mass [m]
length_force = length_force_0;

% Time and initial conditions
time = 0:0.05:10;                                   % Time [s]
x_0 = 0.0;                                          % Initial Condition displacement [m]
x_dot_0 = 0;                                        % Initial Condition velocity [m/s]
phi_0 = 0.0;                                        % Initial Condition angle [rad]
phi_dot_0 = 0;                                      % Initial Condition angle velocity [rad/s]

force                      = 2000;                % Force on the body [N]

%% 2.) Computing
%% 2.) -Numerical solution of the motion
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
        [eigenvector,lambda, ~] = polyeig(C,K,M);
        
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
        x_t = C_vec(1)*exp(lambda(1)*time).*eigenvector(1,1) + C_vec(2)*exp(lambda(2)*time).*eigenvector(1,2) + C_vec(3)*exp(lambda(3)*time).*eigenvector(1,3) +...
            C_vec(4)*exp(lambda(4)*time).*eigenvector(1,4) + x_star(1)*exp(1i*omega*time)+conj(x_star(1))*exp(-1i*omega*time);
        phi_t = C_vec(1)*exp(lambda(1)*time).*eigenvector(2,1) + C_vec(2)*exp(lambda(2)*time).*eigenvector(2,2) + C_vec(3)*exp(lambda(3)*time).*eigenvector(2,3) +...
            C_vec(4)*exp(lambda(4)*time).*eigenvector(2,4)+ x_star(2)*exp(1i*omega*time)+conj(x_star(2))*exp(-1i*omega*time);
        
        % Rewrite homogeneous solution
        x_t = real(x_t);         % Car body displacement
        phi_t = real(phi_t);         % Car body rotation

% Define objective
cost = sum(abs(x_t) + abs(phi_t));  % Define the cost function. Here we will use the sum of the absolute value of the displacement plus the absolute value of phi without any weighing factors.

end    
