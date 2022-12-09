
clearvars; close all; clc;

addpath('../subsystems/');
addpath('../data/');
addpath('../functions/');


configs;

%% Nadir pointing
w0 = [1.5e-3 1.5e-3 1.5e-3]';                % initial angular velocity in body frame [rad/s]

% initial attitude
A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];

q4 = 0.5*sqrt(1 + A_BN0(1, 1) + A_BN0(2, 2) + A_BN0(3, 3));
q1 = 1/(4*q4) * (A_BN0(2, 3) - A_BN0(3, 2));
q2 = 1/(4*q4) * (A_BN0(3, 1) - A_BN0(1, 3));
q3 = 1/(4*q4) * (A_BN0(1, 2) - A_BN0(2, 1));

settings.q0 = [q1 q2 q3 q4]';              % initial estimamted quaternion [-]

n = sqrt(settings.mu/settings.a^3);
A = [                    0            0   ( sat.I(2)-sat.I(3) )/sat.I(1)*n;...
                         0            0                       0;...
     ( sat.I(1)-sat.I(2) )/sat.I(3)*n 0                       0 ];

C = eye(3);
B = diag([1/sat.I(1) 1/sat.I(2) 1/sat.I(3)]);
D = zeros(3);

Kp = [0.9 0.05 0.05]';
Kd = [0.1 0.02 0.02]';