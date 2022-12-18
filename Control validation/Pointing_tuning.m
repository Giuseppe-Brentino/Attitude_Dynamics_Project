
clearvars; close all; clc;

addpath('../subsystems/');
addpath('../data/');
addpath('../functions/');


configs;

%% Nadir pointing
w0 = [1e-2 1e-2 1e-2]';                % initial angular velocity in body frame [rad/s]

% initial attitude
A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];
n = sqrt(settings.mu/settings.a^3);
% [E1,E2,E3] = dcm2angle(A_BN0);
% A0 = angle2dcm(E1+0.7,E2+0.4,E3-1);
A0 = A_BN0;
A = [                    0            0   ( sat.I(2)-sat.I(3) )/sat.I(1)*n;...
                         0            0                       0;...
     ( sat.I(1)-sat.I(2) )/sat.I(3)*n 0                       0 ];

C = eye(3);
B = diag([1/sat.I(1) 1/sat.I(2) 1/sat.I(3)]);
D = zeros(3);

