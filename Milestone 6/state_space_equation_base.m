function dw = state_space_equation_base(time,w,mass,stiffness_f,stiffness_r,damping_f,damping_r,inertia,length_f,length_r,omega,s_max,delta)
%% DESCRIPTION
%
% This is a function for the state space representation of the 2 dof system
% for base excitation
%
%% INPUT
% time:        Time vector
% w:           Initial state vector
% mass:        Mass of car
% stiffness_f: Stiffness front spring
% stiffness_r: Stiffness rear spring
% damping_f:   Stiffness front damper
% damping_r:   Stiffness rear damper
% inertia:     Inertia
% length_f:    Distance front spring damper to center of gravity
% length_r:    Distance rear spring damper to center of gravity
% omega:       Angular frequency of force excitation
% s_max:       Amplitude of base excitation
% delta:       Phase shift between wheels
%
%% OUTPUT
% dw: state space equation
%
%% VERSION
%             author: Martin Lankers (Martin@Lankers.de)
%          copyright: 2017 Lankers Inc.
%      creation date: 03-May-2017
%     Matlab version: 2017a
%            version: 1.0
%
%% REVISION
%
% V1.0 | 03-May-2017 | Martin Lankers | creation
%
%% Compute

% Set up system matrix
A = [0,0,1,0;
    0,0,0,1;
    (-stiffness_r-stiffness_f)/mass,(-stiffness_f*length_f+stiffness_r*length_r)/mass,(-damping_r-damping_f)/mass,(-damping_f*length_f+damping_r*length_r)/mass
    (stiffness_r*length_r-stiffness_f*length_f)/inertia,(-stiffness_f*length_f^2-stiffness_r*length_r^2)/inertia,(-damping_f*length_f+damping_r*length_f)/inertia,(-damping_f*length_f^2-damping_r*length_r^2)/inertia];

% Set up force 1
F1 = stiffness_r*cos(omega*time-delta)-damping_r*omega*sin(omega*time-delta)+stiffness_f*cos(omega*time)-damping_f*omega*sin(omega*time);

% Set up force 2
F2 = (-stiffness_r*cos(omega*time-delta)+damping_r*omega*sin(omega*time-delta))*length_r+(stiffness_f*cos(omega*time)-damping_f*omega*sin(omega*time))*length_f;

% Combine force 1 and force 2 to excitation vector
B = [0;0;(F1*s_max)/mass;(F2*s_max)/inertia];

% Create state space representation
dw = A*w+B;

end