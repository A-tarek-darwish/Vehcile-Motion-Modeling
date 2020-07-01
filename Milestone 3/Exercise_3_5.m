%% DESCRIPTION
%
% This is a Script to solve the differential equation of a force excited single mass 
% system.
%
%% OUTPUT
%
% Formatted figure of the displacement of a single mass system and its
% animation.
%
%
%% Program
clear                                             % Delete Workspace
clc                                               % Clear Command Window
close all                                         % Close all figures

%% 1.) Definitions
%% 1.) -Parameter definition
mass                       = 750;                 % Mass of the body [kg]
stiffness                  = 50000;               % Stiffness Coefficient of spring [N/m]
damping                    = 10;                   % Damping coefficient of damper [Ns/m]
time                       = 0:0.01:10;            % Time [s]

x_0                        = 0.00;                % Initial Condition displacement
x_dot_0                    = 0.00;                   % Initial Condition velocity

mass_extruded              = 0.03;                % Mass of the excitation [kg]
omega                      = 6.28;                % Angular frequency of the excitation [1/s]
radius                     = 0.24;                % Radius of the excitation [m]

force = omega^2*radius*mass_extruded;             % Calculate the force with given parameters

n = 10;                                           % Number of test calculations
    
%% 2.) Computing
%%%%%%%% analytically %%%%%%%%
for k = 1:n
    tic
    %% 2.) -analytically | Parameter calculation
    dampingcoefficient = damping/(2*mass);              % Calculate damping coefficient
    angulareigenfrequency = sqrt(stiffness/mass);       % Calculate angular eigenfrequency
    
    %% 2.) -analytically | Calculation of the characteristic polynomial
    lambda = roots([1, 2*dampingcoefficient, angulareigenfrequency^2]); %Calculate roots
    
    %% 2.) -analytically | Calculation of the constants
    x_max = sqrt(force^2/((damping*omega)^2 + (stiffness-mass*omega^2)^2));     % maximum amplitude of particular solution
    phi = atan((damping*omega)/(stiffness-mass*omega^2));                  % phase shift of the particular solution
    
    k1 = x_0 - (lambda(1)*x_0 - x_dot_0 - x_max*omega*sin(-phi) - x_max*lambda(1)*cos(-phi)) / (lambda(1) - lambda(2))- x_max*cos(-phi);                                            % Calculation of the constant 1
    k2 = (lambda(1)*x_0 - x_dot_0 - x_max*omega*sin(-phi) - x_max*lambda(1)*cos(-phi)) / (lambda(1) - lambda(2)); % Calculation of the constant 2
    
    %% 2.) -analytically | Calculation of the solution
    x_t_h = k1 * exp(lambda(1)*time) + k2 * exp(lambda(2)*time);                         % Evaluate function at points "time"
    v_t_h = k1 * lambda(1) * exp(lambda(1)*time) + k2 * lambda(2) * exp(lambda(2)*time); % Evaluate function at points "time"
    
    x_t_p = x_max*cos(omega*time-phi);                                          % x(t) of the particular solution
    v_t_p = -x_max*omega*sin(omega*time-phi);                                   % v(t) of the particular solution
    
    x_t = real(x_t_h) + x_t_p;            % Calculate overall solution by superposition of homogeneous and particular solution
    v_t = real(v_t_h) + v_t_p;            % Calculate overall solution by superposition of homogeneous and particular solution
    time_analytically(k) = toc;
end

%%%%%%%% dsolve %%%%%%%%
%% 2.) -dsolve | Calculation of the solution
for k = 1:n
    tic
    t = sym('t');                                     % Define time variable
    x = symfun(sym(str2sym('x(t)')), [t]);            % Define dependent variable
    Dx = diff(x,1);                                   % Define first derivation
    D2x = diff(x,2);                                  % Define second derivation
    x = dsolve( mass*D2x + damping*Dx + stiffness*x == force*cos(omega*t)  ,x(0)==x_0  ,Dx(0) == x_dot_0,  't');  % Initial conditions units are [m] and [m/s] respectively
        
    %% 2.) -dsolve | Evaluate the equation
    x_fun = matlabFunction(x);                          % Create function handle for x
    x_dot_fun = matlabFunction(diff(x));                % Create function handle for x_dot    
    time_dsolve(k) = toc;

end

%%%%%%%% ode %%%%%%%%
%% 2.) -ode | Set initial conditions
for k = 1:n
    tic
    w0 = [x_0 , x_dot_0];                             % Create a vector with initial conditions
    
    %% 2.) -ode | Numerical solution of the motion
    A = [0, 1;  (-1)*stiffness / mass, (-1)*damping / mass ]; % Create system Matrix
    dw = @(t,w) A*w;    % Define derivative
    [tsim,wsim] = ode45(dw,time,w0); % Calling numerical solver
    time_ode(k) = toc;
        
end
%% 3.) Plot
compare_analytical_dsolve = sum(time_dsolve)/sum(time_analytically);
compare_analytical_ode = sum(time_ode)/sum(time_analytically);
compare_dsolve_ode = sum(time_dsolve)/sum(time_ode);


