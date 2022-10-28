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

% initial attitude and position conditions (Ephemeris 00:00:00 27-10-2022)

settings.a = 7099.67959313555;              % major semi-axis [km] 
settings.e = 0.00117256328941186;           % eccentricity [-]
settings.i = deg2rad(91.9344391020783);     % inclination [rad]
settings.OM = deg2rad(315.709234758653);    % RAAN [rad]
settings.om = deg2rad(101.432935891144);    % pericenter's anomaly [rad]
settings.theta = 0;                         % true anomaly [rad] (assumed equal to zero)

settings.mu = astroConstants(13);            % Earth's planetary constant [km^3/s^2]

[r0, v0] = kep2car(settings);    % initial position and velocity in inertial frame [km]
settings.w0 = [];               % initial angular velocity in body frame [rad/s]
settings.E0 = [];                % initial euler angles [rad]

% environment parameters

environment.E_w =  deg2rad(15.04/3600);         % Earth's angular velocity around it's axis [rad/s]
environment.theta_G0 = 0;                       % Initial Greenwich Meridian Longitude [rad] (DA CAMBIARE)
environment.theta_E0 = 0;                       % Initial Earth longitude wrt the sun [rad]  (DA CAMBIARE)
environment.S_w = 2*pi/(365.25*24*3600);        % Earth's angular velocity around the sun [rad/s]
environment.date = 221*24*3600;                 % Seconds from spring equinox (March 20) to Oct 27 [s]
environment.Rs = astroConstants(2);             % Sun-Earth distance (1 AU) [km]   
environment.eps = deg2rad(23.45);               % Inclination of the ecliptic plane [rad]

% satellite parameters
sat.Ix = [];
sat.Iy = [];
sat.Iz = [];

% sensor parameters
sensor.boh = [];

% actuator parameters
coils.boh = [];
secondActuator.boh = [];

% test = sim("orbit_propagation.slx");


