%% Main script of the Spacecraft Attitude Dynamics course

%{
    Specifications:
    - Large satellite ( over 500kg )
    - Euler angles
    - Star sensor
    - 3 Magnetic coils 

%}

clearvars; close all; clc;

addpath('./subsystems/');
addpath('./data/');
addpath('./functions/');


configs;

settings.w0 = sqrt(settings.mu/settings.a^3)*[0 1 0]';     % initial angular velocity in body frame [rad/s]

% initial attitude
% A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];
h0=cross(r0,v0);
A_BN0 = [cross(-r0',h0')/norm(cross(r0,h0)); h0'/norm(h0); r0'/norm(r0)];

% Euler angles DA CONTROLLARE!!!
% phi0=atan(-A_BN0(2,1)/A_BN0(2,2)); 
% theta0=asin( A_BN0(2,3)); 
% psi0=atan(-A_BN0(1,3)/A_BN0(3,3)); 

q4 = 0.5*sqrt(1 + A_BN0(1, 1) + A_BN0(2, 2) + A_BN0(3, 3));
q1 = 1/(4*q4) * (A_BN0(2, 3) - A_BN0(3, 2));
q2 = 1/(4*q4) * (A_BN0(3, 1) - A_BN0(1, 3));
q3 = 1/(4*q4) * (A_BN0(1, 2) - A_BN0(2, 1));

[E1, E2, E3] = dcm2angle(A_BN0, "ZXY");

settings.E0 = [E1 E2 E3]';                % initial euler angles [rad]
settings.q0 = [q1 q2 q3 q4]';             % initial estimamted quaternion [-]
