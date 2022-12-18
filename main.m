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
addpath('./Control validation/');


configs;

settings.w0 = [1e-2, 1e-2, 1e-2]';  % detumbling initial angular velocity in body frame [rad/s] = 3 deg/s
% settings.w0 = [0.035 0.035 0.035]';     % initial angular velocity in body frame [rad/s] = 2 deg/s
% settings.w0 = [0.005 0.005 0.005]';     % initial angular velocity in body frame [rad/s] = 0.3 deg/s
% settings.w0 = n;                        % star tracker test angular velocity
% settings.w0 = [deg2rad(1), n(2), 0];

% initial attitude
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
