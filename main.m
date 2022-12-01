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
% initial attitude and position conditions (Ephemeris 00:00:00 27-10-2022)

settings.a = 7099.67959313555;              % major semi-axis [km] 
settings.e = 0.00117256328941186;           % eccentricity [-]
settings.i = deg2rad(91.9344391020783);     % inclination [rad]
settings.OM = deg2rad(315.709234758653);    % RAAN [rad]
settings.om = deg2rad(101.432935891144);    % pericenter's anomaly [rad]
settings.theta = 0;                         % true anomaly [rad] (assumed equal to zero)

settings.mu = astroConstants(13);           % Earth's planetary constant [km^3/s^2]

[r0, v0] = kep2car(settings);    % initial position and velocity in inertial frame [km]
settings.w0 = [0 0 0]';                % initial angular velocity in body frame [rad/s]
settings.E0 = [0 0 0]';                % initial euler angles [rad]
settings.q0 = [0 0 0 1]';              % initial estimamted quaternion [-]

settings.N = 13;                 % order of the magnetic field model [-]
% environment parameters

environment.E_w =  deg2rad(15.04/3600);         % Earth's angular velocity around it's axis [rad/s]
environment.theta_G0 = 0;                       % Initial Greenwich Meridian Longitude [rad] (mezzanotte 21Settembre)
environment.theta_E0 = 0;                       % Initial Earth longitude wrt the sun [rad]  (DA CAMBIARE)
environment.S_w = 2*pi/(365.25*24*3600);        % Earth's angular velocity around the sun [rad/s]
environment.date = 221*24*3600;                 % Seconds from spring equinox (March 20) to Oct 27 [s]
environment.Rs = astroConstants(2);             % Sun-Earth distance (1 AU) [km]   
environment.eps = deg2rad(23.45);               % Inclination of the ecliptic plane [rad]
environment.R = 6371.2;                         % Earth's radius [km]

load("WMM.mat");                                % World Magnetic Model parameters 
environment.WMM.K = K;                          % [-]
environment.WMM.g = g;                          % [nT]
environment.WMM.h = h;                          % [nT]

% satellite parameters
sat.m  = 720;                                   % Spacecraft mass [kg]
sat.dipole = 3.5e-3*sat.m * ones(3,1);          % Spacecraft dipole moment [Am^2]
sat.I = [1000 500 1500]';                               % Column vector with Principal Inertia Moments
theta = deg2rad(38.35);                         % angle between z_body and solar panel DA TROVARE SUL CAD[rad]
sat.panel1 = [0;-sin(theta);cos(theta)];        % normal to 1st solar panel in body frame
sat.panel2 = [0;sin(theta);cos(theta)];         % normal to 2nd solar panel in body frame

% sensor parameters
load("star_catalogue.mat")                     
sensors.star.fov = 20;                           % nominal field of view [deg]
sensors.star.alpha_err = (deg2rad(2/3600))^2;    % variance of alpha [rad]
sensors.star.delta_err = (deg2rad(4.91/3600))^2; % variance of delta [rad]
sensors.star.frequency = 0.3;                      % maximum update rate [Hz]
sensors.star.inclination = deg2rad(15);          % inclination of the sensor wrt z_body axis [rad]
sensors.star.focal_length = 20;                  % focal length of the sensor [mm]
sensors.star.pixel = sensors.star.focal_length...% pixel length [mm]
    *tan(deg2rad(sensors.star.fov/2))/512;
theta = pi/2 - sensors.star.inclination;
sensors.star.S = zeros(4, 10);  
sensors.star.alpha = zeros (10, 1); 
sensors.star.delta = zeros (10, 1); 
sensors.star.ASB1 = [ 1       0          0; ...         % rotation matrix body to sensor 1
                      0  cos(theta)  sin(theta); ...
                      0  -sin(theta) cos(theta)];
sensors.star.ASB2 = [ 1       0          0; ...
                      0  cos(-theta)  sin(-theta); ...  % rotation matrix body to sensor 2 CHECK
                      0  -sin(-theta)  cos(-theta)];
sensors.star.opt_w = deg2rad(0.3);
sensors.star.reduced_perf_w = deg2rad(1.5);
% magnetic sensor parameters
sensors.mag.acc = 0.5*1e-2;                              % accuracy percentage of full scale
sensors.mag.lin = 0.015*1e-2;                            % linearity percentage of full scale
sensors.mag.FMR = 100*1e3;                                 % field measurement of range [nT]
sensors.mag.STD_dev = sqrt( (sensors.mag.acc/(sqrt(3)))^2 ...
    + (sensors.mag.lin/(sqrt(3)))^2 )*sensors.mag.FMR;        % standard deviation 
sensors.mag.variance = (sensors.mag.STD_dev)^2;            % variance 
sensors.mag.sensitivity = 100*1e-6;                        % sensitivity [V/nT]
sensors.mag.quant_V = 0.050;                                 % quantization interval [V]
sensors.mag.quant_T = sensors.mag.quant_V/sensors.mag.sensitivity;      % quantization interval [nT]
sensors.mag.frequency = 1/5;                                    % sample time [s]
sensors.mag.angles.ax=deg2rad(0.5);
sensors.mag.angles.ay=deg2rad(0.5);
sensors.mag.angles.az=deg2rad(0.5);
sensors.mag.angles.bx=deg2rad(10);
sensors.mag.angles.by=deg2rad(20);
sensors.mag.angles.bz=deg2rad(30);

ax=sensors.mag.angles.ax;
ay=sensors.mag.angles.ay;
az=sensors.mag.angles.az;
bx=sensors.mag.angles.bx;
by=sensors.mag.angles.by;
bz=sensors.mag.angles.bz;

A_mag=[cos(ax) sin(ax)*cos(bx) sin(ax)*sin(bx);
    sin(ay)*cos(by) cos(ay) sin(ay)*sin(by);
    sin(az)*cos(bz) sin(az)*sin(bz) cos(az)];

sensors.mag.A_mag = A_mag;

% actuator parameters

thruster.thrust = 0.01;                                 % Nominal thrust [N]   
thruster.direction = [ 0  0  0  0  0  0  0  0; ...      % Direction in which the thrust is applied
                       0  0  0  0  1 -1 -1  1; ...
                       1 -1  1 -1  0  0  0  0; ]; 
thruster.position = ones(3,8);                          % Thruster position wrt center of mass [m]
thruster.firing_time = 0.05;                            % Minimum firing time [s]

magnetorquers.dipole = 30;                              % Maximum magnetic dipole [Am^2] da datsheet (attuatore su eoportal)

% test = sim("orbit_propagation.slx");


