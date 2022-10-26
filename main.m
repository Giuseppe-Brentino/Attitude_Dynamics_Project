%% Main script of the Spacecraft Attitude Dynamics course

%{
    Specifications:
    - Large satellite ( over 500kg )
    - Euler angles
    - Star sensor
    - 3 Magnetic coils 

%}

clearvars; close all; clc;

% initial attitude and position conditions
r0 = [];                        % initial position in inertial frame [km]
v0 = [];                        % initial velocity in inertial frame [km/s]

w0 = [];                        % initial angular velocity in body frame [rad/s]
E0 = [];                        % initial euler angles [rad]

% satellite parameters
sat.Ix = [];
sat.Iy = [];
sat.Iz = [];

% sensor parameters
sensor.boh = [];

% actuator parameters
coils.boh = [];
secondActuator.boh = [];


