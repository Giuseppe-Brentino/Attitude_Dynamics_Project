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

create_plots = input('Do you want to generate the plots? ("yes" or "no"): ','s');

if strcmp(create_plots,'yes')
    save_plots = input('Do you want to save the plots in pdf format? ("yes" or "no"): ','s');
end

configs;

%% Initial conditions after release from launcher

% initial angular velocity in body frame [rad/s]
settings.w0 = [3e-2,1e-2,3e-2]'; 

% initial attitude
A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];

% initial Euler angles
phi0= pi + atan(-A_BN0(2,1)/A_BN0(2,2));
theta0=asin( A_BN0(2,3)); 
psi0=atan(-A_BN0(1,3)/A_BN0(3,3)); 
settings.E0 = [phi0,theta0,psi0]';                % initial euler angles [rad]
clear phi0 theta0 psi0

% initial quaternion guess
q4 = 0.5*sqrt(1 + A_BN0(1, 1) + A_BN0(2, 2) + A_BN0(3, 3));
q1 = 1/(4*q4) * (A_BN0(2, 3) - A_BN0(3, 2));
q2 = 1/(4*q4) * (A_BN0(3, 1) - A_BN0(1, 3));
q3 = 1/(4*q4) * (A_BN0(1, 2) - A_BN0(2, 1));
settings.q0 = [q1 q2 q3 q4]';             % initial estimamted quaternion [-]
clear q1 q2 q3 q4

%% De-tumbling

control.algorithm = 'De-tumbling';
magnetorquers.dipole = 120;               % Maximum magnetic dipole [Am^2]
phase1 = sim('Model.slx');
save phase1 phase1
%% Pointing

% initial conditions ( final conditions of phase 1 )
environment.theta_G0 = phase1.theta_G(end); 
r0 = phase1.r_N(end,:)';
v0 = phase1.v(end,:)';

settings.w0 = deg2rad(phase1.w_BN(end,:))';
settings.E0 = deg2rad(phase1.E_312(end,:)');

control.algorithm = 'Pointing';
magnetorquers.dipole = 350;      % Maximum magnetic dipole [Am^2]
phase2 = sim('Model.slx');

save phase2 phase2

%% plots
if strcmp(create_plots,'yes')
    plots;
end
