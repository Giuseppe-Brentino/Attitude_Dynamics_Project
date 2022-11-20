%% Main script of the Spacecraft Attitude Dynamics course
% Modified for testing results

% NOTE: loro usano come initial condition A_0=eye(3), noi useremo Euler
%       Angles E_0=[0,0,0]';

clearvars; close all; clc;

addpath('./subsystems/');
addpath('./data/');
addpath('./functions/');

%% Parametri copiati dallo script del quartetto maggico
G = 6.67408e-20;
mt = 5.97219e24;
r_gc = 0.015*[1,1,1];  % Che roba è???

orbit = 'GEO';

switch orbit
    case 'LEO'
        a = 6373+400;
        e = 0.3;
        i = deg2rad(20);
        n=sqrt(G*mt/a^3);
        SimTime = 2*(2*pi/n);
    case 'GEO'
        a = 6373+35785;
        e = 0;
        i = deg2rad(0);
        n=sqrt(G*mt/a^3);
        SimTime = 2*pi/n;
    case 'other'
        a = 2.6e4;
        e = 0.0618;
        i = deg2rad(0);
        n = sqrt(G*m/a^3);
        SimTime = 2*pi/n;
end

m = [0.01, 0.05, 0.01]'; %[A m2] Parasitic moment

A_0 = eye(3);

I_x = 1.009;
I_y = 0.251;
I_z = 0.916;


%% initial attitude and position conditions (Ephemeris 00:00:00 27-10-2022)
settings.a = a;                                                             % major semi-axis [km] 
settings.e = e;                                                             % eccentricity [-]
settings.i = i;                                                             % inclination [rad]
% Other parameters supposed by Nic Galletta
settings.OM = deg2rad(0);                                                   % RAAN [rad]
settings.om = deg2rad(0);                                                   % pericenter's anomaly [rad]
settings.theta = 0;                                                         % true anomaly [rad] (assumed equal to zero)
 
settings.mu = astroConstants(13);                                           % Earth's planetary constant [km^3/s^2]

%% Initial conditions 
[r0, v0] = kep2car(settings);                                               % initial position and velocity in inertial frame [km]
settings.w0 = [0 0 n]';                                                     % initial angular velocity in body frame [rad/s]
settings.E0 = [0 0 0]';                                                     % initial euler angles [rad]

settings.N = 13;                                                            % order of the magnetic field model [-]

%% environment parameters
environment.E_w =  deg2rad(15.04/3600);                                     % Earth's angular velocity around it's axis [rad/s]
environment.theta_G0 = 0;                                                   % Initial Greenwich Meridian Longitude [rad] (mezzanotte 21Settembre)
environment.theta_E0 = 0;                                                   % Initial Earth longitude wrt the sun [rad]  (DA CAMBIARE)
environment.S_w = 2*pi/(365.25*24*3600);                                    % Earth's angular velocity around the sun [rad/s]
environment.date = 221*24*3600;                                             % Seconds from spring equinox (March 20) to Oct 27 [s]
environment.Rs = astroConstants(2);                                         % Sun-Earth distance (1 AU) [km]   
environment.eps = deg2rad(23.45);                                           % Inclination of the ecliptic plane [rad]
environment.R = 6371.2;                                                     % Earth's radius [km]

load("WMM.mat");                                                            % World Magnetic Model parameters 
environment.WMM.K = K;                                                      % [-]
environment.WMM.g = g;                                                      % [nT]
environment.WMM.h = h;                                                      % [nT]

%% satellite parameters
%sat.m  = 720;                                                               % Spacecraft mass [kg]
%sat.dipole = 3.5e-3*sat.m * ones(3,1);                                      % Spacecraft dipole moment [Am^2]
sat.dipole = m;                                                             % Spacecraft dipole moment [Am^2]
sat.I = [I_x I_y I_z]';                                                     % Column vector with Principal Inertia Moments
theta = deg2rad(70);                                                        % angle between z_body and solar panel DA TROVARE SUL CAD[rad]
sat.panel1 = [0;-sin(theta);cos(theta)];                                    % normal to 1st solar panel in body frame
sat.panel2 = [0;sin(theta);cos(theta)];                                     % normal to 2nd solar panel in body frame

% sensor parameters
load("star_catalogue.mat")
sensors.star.Bias_max = 10;                                                 % max bias error [arcsec] 
sensors.star.fov = 20;                                                      % nominal field of view [deg]
sensors.star.spatial_err = deg2rad(1.5/3600);                               % star position error in the sensor at 3 sigma [arcsec]
sensors.star.frequency = 5;                                                 % maximum update rate [Hz]
sensors.star.inclination = deg2rad(15);                                     % inclination of the sensor wrt z_body axis [rad]
sensors.star.focal_length = 20;                                             % focal length of the sensor [mm]
sensors.star.pixel = sensors.star.focal_length...                           % pixel length [mm]
    *tan(deg2rad(sensors.star.fov/2))/512;
theta = pi/2 - sensors.star.inclination;
sensors.star.S = zeros(4, 10);  
sensors.star.alpha = zeros (10, 1); 
sensors.star.delta = zeros (10, 1); 
sensors.star.ASB1 = [ 1       0          0; ...                             % rotation matrix body to sensor 1
                      0  cos(theta)  sin(theta); ...
                      0  -sin(theta) cos(theta)];
sensors.star.ASB2 = [ 1       0          0; ...
                      0  cos(-theta)  sin(-theta); ...                      % rotation matrix body to sensor 2
                      0  sin(-theta)  -cos(-theta)];

% magnetic sensor parameters
sensors.mag.acc = 0.5*1e-2;                                                 % accuracy percentage of full scale
sensors.mag.lin = 0.015*1e-2;                                               % linearity percentage of full scale
sensors.mag.FMR = 100*1e3;                                                  % field measurement of range [nT]
sensors.mag.STD_dev = sqrt( (sensors.mag.acc/(sqrt(3)))^2 ...
    + (sensors.mag.lin/(sqrt(3)))^2 )*sensors.mag.FMR;                      % standard deviation 
sensors.mag.variance = (sensors.mag.STD_dev)^2;                             % variance 
sensors.mag.sensitivity = 100*1e-6;                                         % sensitivity [V/nT]
sensors.mag.quant_V = 0.050;                                                % quantization interval [V]
sensors.mag.quant_T = sensors.mag.quant_V/sensors.mag.sensitivity;          % quantization interval [nT]
sensors.mag.frequency = 1/5;                                                % sample time [s]

% actuator parameters
coils.boh = [];
secondActuator.boh = [];

test = sim("Model.slx");


%% Plot
tvec=test.tout;

figure
plot(test.w_BN)
xlim([tvec(1) tvec(end)]); grid on;
title('w_BN')

figure
plot(test.A_BN)
xlim([tvec(1) tvec(end)]); grid on;
title('A_BN')

figure
plot(test.A_BL)
xlim([tvec(1) tvec(end)]); grid on;
title('A_BL')

figure
plot(test.A_LN)
xlim([tvec(1) tvec(end)]); grid on;
title('A_LN')

figure
plot(test.GG_Torque)
xlim([tvec(1) tvec(end)]); grid on;
title('Gravity gradient disturbance torque')

figure
plot(test.Magnetic_torque)
xlim([tvec(1) tvec(end)]); grid on;
title('Magnetic disturbance torque')

figure
plot(test.tot_torque)
xlim([tvec(1) tvec(end)]); grid on;
title('Total torque')
