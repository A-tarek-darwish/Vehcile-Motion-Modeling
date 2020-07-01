%% DESCRIPTION
%
% This is a Script to animate a two degree of freedom system under base excitation.

%
%% OUTPUT
%
% Formatted figure of the displacement and angle of a two degree of freedom system and its
% animation.
%% 1.) Definitions
%% 1.) -Parameter definition
% Masses and inertias
mass                      = 1000;                   % Mass of the body [kg]
inertia                   = 500;                   % Inertia of the body [kg*m^2]

% Stiffness and damping values
stiffness_f               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_f                 = 500;                   % Damping coefficient of damper [Ns/m]
stiffness_r               = 60000;                  % Stiffness coefficient of spring [N/m]
damping_r                 = 500;                   % Damping coefficient of damper [Ns/m]

% Lengths center of gravity to front and rear end
length_f                  = 2.5;                    % Distance of the right spring-damper to the center of mass [m]
length_r                  = 2.5;                    % Distance of the left spring-damper to the center of mass [m]

% Time and initial conditions
tend = 15;
time = 0:0.01:tend;                                   % Time [s]

x_0 = 0.0;                                          % Initial Condition displacement [m]
x_dot_0 = 0;                                        % Initial Condition velocity [m/s]
phi_0 = 0.0/180*pi;                                 % Initial Condition angle [rad]
phi_dot_0 = 0;                                      % Initial Condition angle velocity [rad/s]
w0 = [x_0,phi_0,x_dot_0,phi_dot_0];                               % Create a vector with initial conditions

velocity = 50/3.6;             % Velocity of car in km/h needs to be converted into m/s
s_max = 0.01;                  % Amplitude of base excitation is 10 mmm

length_exc_lower_bound = 2;                       % Define a lower bound of the wavelength of the excitation of 2 
length_exc_upper_bound = 50;                       % Define a upper bound of the wavelength of the excitation of 50
length_exc_vector = length_exc_lower_bound:0.25:length_exc_upper_bound; % Define the length_exc_vector from upper bound to lower bound with an inkrement of 0.25
k = 0;          % Initialise counting variable
for length_exc = length_exc_vector  % Loop over length_exc_vector
    k = k + 1;  % Counting variable goes one up
    omega(k) = (velocity*2*pi)/length_exc;                              % Angular Velocity of excitation     
    
    delta = ((length_f+length_r)*2*pi)/length_exc;                   % Phase shift between front and rear wheel
    phi_0 = asin(abs((s_max*cos(-delta)-s_max))/(length_f+length_r))+phi_0;     % Initial consistency condition to avoid step excitation (rear wheel)
    x_0 = -length_r*abs((s_max*cos(-delta)-s_max))/(length_f+length_r)+ s_max +x_0;   % Initial consistency condition to avoid step excitation (rear wheel)
    
    %% 2.) Computing
    %% 2.) -Numerical solution of the motion
    [tsim,wsim{k}] = ode45(@state_space_equation_base,time,w0,'options',mass,stiffness_f,stiffness_r,damping_f,damping_r,inertia,length_f,length_r,omega(k),s_max,delta); % Calling numerical solver
    
end

excitation_frequency = omega/2/pi;  % Calculate the excitation frequency from the angular velocity of excitation, remember omega = 2*pi*excitation_frequency

% Now we want to plot some of the results in the time domain to have a look
% at the displacement 
u = 0;   % Initialize counting variable
figure 
for n = [1 3 4 7]  % Put in some numbers of simulation results 
 plot(tsim,wsim{n}(:,1))  % Plot the results of the displacement against the time vector
 hold on  % You want every result in just one plot 
 u = u +1; % The counting variable goes one up
 lgnd{u} = sprintf('F %.2f Hz',excitation_frequency(u));  % You want to generate a legend entry for every loop which gives you the information about the according excitation frequency
end
legend(lgnd,'Location','northeastoutside')  % Plot the legend in the location northeastoutside

% As you have seen in the time result plots, there is a startup tranisient response, but we want only the amplitude of the steady state response (particular solution).
% Therefore we just want to analyze the last 10 percent of the signal.
number_of_timesteps = length(tsim);  % Determine the number of total timesteps first.
number_of_timesteps_to_analyze = floor(0.1*number_of_timesteps);   % Define the number of timesteps that you want to analyze. We want to analyze only 10 percent of the total time steps. If you get a non-integer number, round off to the next integer.

idx_timesteps_to_analyze = (number_of_timesteps-number_of_timesteps_to_analyze):number_of_timesteps;  % Determine the indices of the timesteps to analyze 
x_hat = cellfun(@(x) max(abs(x(idx_timesteps_to_analyze,1))),wsim);  % Analyze every simulation and find the maximum of the absolute values of the displacement. You could do that either with a for loop or with a cell function. type help cellfun into the command window  

figure
plot(excitation_frequency,x_hat)  % Plot the amplitude over the excitation frequency
xlabel('[Hz]')                  % Label the x axis with '[Hz]'
ylabel('[m]')           % Label the y axis with '[m]'

[value_max,idx_max] = max(x_hat);  % Find the maximum value of all amplitudes and its idx
length_exc_critical = length_exc_vector(idx_max);  % Find the worst road conditions by determining the critical wavelength
