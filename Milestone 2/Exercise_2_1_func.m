%% VERSION
%             author: Martin Lankers (Martin.Lankers.de)
%      creation date: 25-April-2018
%     Matlab version: 2017a
%
%% REVISION
%
% V1.0 | 25-April-2018 | Martin Lankers | creation

function [x_t,v_t] = Exercise_2_1_func(mass,stiffness,damping,time,x_0,x_dot_0)
    %% 2.) Computing
    %% 2.) -Parameter calculation
    dampingcoefficient = damping/(2*mass);              % Calculate damping coefficient
    angulareigenfrequency = sqrt(stiffness/mass);       % Calculate angular eigenfrequency

    %% 2.) -Calculation of the characteristic polynomial
    lambda = roots([1, 2*dampingcoefficient, angulareigenfrequency^2]); %Calculate roots

    %% 2.) -Calculation of the constants
    k1 = (x_dot_0 - lambda(2)*x_0) / (lambda(1) - lambda(2));           % Calculation of the constant 1
    k2 = (lambda(1)*x_0 - x_dot_0) / (lambda(1) - lambda(2));           % Calculation of the constant 2

    %% 2.) -Calculation of the solution
    x_t_h = k1 * exp(lambda(1)*time) + k2 * exp(lambda(2)*time);                         % Evaluate function at points "time"
    v_t_h = k1 * lambda(1) * exp(lambda(1)*time) + k2 * lambda(2) * exp(lambda(2)*time); % Evaluate function at points "time"

    x_t = real(x_t_h);                                                  % Make sure only real part is considered
    v_t = real(v_t_h);                                                  % Make sure only real part is considered
end

