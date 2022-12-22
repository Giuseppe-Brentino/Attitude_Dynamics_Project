
clearvars; close all; clc;

addpath('../subsystems/');
addpath('../data/');
addpath('../functions/');


configs;

%% Nadir pointing
settings.w0 = [1e-2 1e-2 1e-2]';                % initial angular velocity in body frame [rad/s]

% initial attitude
A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];
n = sqrt(settings.mu/settings.a^3);
[E1,E2,E3] = dcm2angle(A_BN0,'ZXY');
A0 = angle2dcm(E1-0.5,E2+0.8,E3+0.3,'ZXY');

% K_w = 1e-2;
% K_cz = [1 1 1]; %sat.I/sat.I(3)*1e-4;
% K_cx = [1 1 1]; %sat.I/sat.I(3)*1e-4;

thruster.K_e = 47;
thruster.K_m = 40;
thruster.T_m = 20;
thruster.Uon = 19;
thruster.Uoff = 18.9;

q4 = 0.5*sqrt(1 + A_BN0(1, 1) + A_BN0(2, 2) + A_BN0(3, 3));
q1 = 1/(4*q4) * (A_BN0(2, 3) - A_BN0(3, 2));
q2 = 1/(4*q4) * (A_BN0(3, 1) - A_BN0(1, 3));
q3 = 1/(4*q4) * (A_BN0(1, 2) - A_BN0(2, 1));
settings.q0 = [q1 q2 q3 q4]';             % initial estimamted quaternion [-]
settings.E0 = [E1-0.5,E2+0.8,E3+0.3]';                % initial euler angles [rad]