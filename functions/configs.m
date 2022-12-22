% initial attitude and position conditions (Ephemeris 00:00:00 27-10-2022)
rng('default');
settings.a = 7099.67959313555;              % major semi-axis [km] 
settings.e = 0.00117256328941186;           % eccentricity [-]
settings.i = deg2rad(91.9344391020783);     % inclination [rad]
settings.OM = deg2rad(315.709234758653);    % RAAN [rad]
settings.om = deg2rad(101.432935891144);    % pericenter's anomaly [rad]
settings.theta = 0;                         % true anomaly [rad] (assumed equal to zero)

settings.mu = astroConstants(13);                       % Earth's planetary constant [km^3/s^2]
settings.Time = 2*pi*sqrt(settings.a^3/settings.mu);    % Orbital period [s]
n = sqrt(settings.mu/settings.a^3);

[r0, v0] = kep2car(settings);    % initial position and velocity in inertial frame [km]


settings.N = 13;                 % order of the magnetic field model [-]

%% environment parameters
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

environment.Fe=1358;                            % Solar Radiation 1AU [W/m2]
environment.c= 299792458;                       % speed of light [m/s]
environment.P=environment.Fe/environment.c;     % Solar pressure [N/m^2]

%% satellite parameters
sat.m  = 720;                                   % Spacecraft mass [kg]
sat.dipole = 3.5e-3*sat.m * ones(3,1);          % Spacecraft dipole moment [Am^2]
sat.I = [279 945 1085]';                        % Column vector with Principal Inertia Moments

% reflectance coefficients
sat.rho_s.panel=0.1;
sat.rho_s.side=0.5;
sat.rho_d.panel=0.1;
sat.rho_d.side=0.1;
alfa = deg2rad(38.35);                          % angle between z_body and normal to solar panel [rad]
beta = deg2rad(8.84);                           % angle between y_body and S/C sides [rad]
sat.Normals = [-1     0          0       0     0       1     0;  ...      % normal versors to each surface in body frame
                0 -sin(alfa) -cos(beta)  0  sin(alfa)  0   cos(beta); ...
                0  cos(alfa) -sin(beta) -1  cos(alfa)  0  -sin(beta)];
sat.rGC = 1e-3* [-1715    88    90   88   88  1802    90; ...   % vectors CoG to centroid of surfaces [m]
                     0  -611 -1113    0  611     0  1113; ...
                   152   618  -119 -419  618   152  -119];
sat.Aree = 1e-6* [2402565 5048400 2187269 7981554 5048400 2402565 2187269];  % area of S/C surfaces [m^2]
sat.rhos_panel=0.5;
sat.rhod_panel=0.1;
sat.rhos=0.1;
sat.rhod=0.1;
sat.panel1=sat.Normals(:,2);
sat.panel2=sat.Normals(:,5);

%% sensor parameters
load("star_catalogue.mat")                     
sensors.star.fov = 20;                           % nominal field of view [deg]
sensors.star.alpha_std = deg2rad(2/3600);    % standard deviation of alpha [rad]
sensors.star.delta_std = deg2rad(4.91/3600); % standard deviation of delta [rad]
sensors.star.alpha_err = sensors.star.alpha_std^2;    % variance of alpha [rad]
sensors.star.delta_err = sensors.star.delta_std^2; % variance of delta [rad]
sensors.star.frequency = 1;                      % maximum update rate [Hz]
% sensors.star.inclination = deg2rad(68.7);          % inclination of the normal of the sensor wrt z_body axis [rad]
sensors.star.inclination = deg2rad(75);          % inclination of the normal of the sensor wrt z_body axis [rad]
sensors.star.focal_length = 20;                  % focal length of the sensor [mm]
sensors.star.pixel = sensors.star.focal_length...% pixel length [mm]
    *tan(deg2rad(sensors.star.fov/2))/512;
sensors.star.S = zeros(4, 10);  
sensors.star.alpha = zeros (10, 1); 
sensors.star.delta = zeros (10, 1); 
theta = sensors.star.inclination;
sensors.star.ASB1 = [ 1       0          0; ...         % rotation matrix body to sensor 1
                      0  cos(theta)  sin(theta); ...
                      0  -sin(theta) cos(theta)];
sensors.star.ASB2 = [ 1       0          0; ...
                      0  cos(-theta)  sin(-theta); ...  % rotation matrix body to sensor 2 CHECK
                      0  -sin(-theta)  cos(-theta)];
sensors.star.opt_w = deg2rad(0.3);
sensors.star.reduced_perf_w = deg2rad(1.5);
alfa1 = sensors.star.inclination - alfa - deg2rad(sensors.star.fov/2);   % lower angle between solar panel normal and star tracker normal 
                                                                         % such that sun is in fov of star tracker 
alfa2 = sensors.star.inclination - alfa + deg2rad(sensors.star.fov/2);   % greater angle between solar panel normal and star tracker normal
                                                                         % such that sun is in fov of star tracker 
sensors.star.startup = 1800;                        % seconds after which it's possible to use the star sensor

% magnetic sensor parameters
sensors.mag.acc = 0.5;                              % accuracy percentage of full scale
sensors.mag.lin = 0.015;                            % linearity percentage of full scale
sensors.mag.FMR = 100*1e-6;                                 % field measurement of range [T]
sensors.mag.STD_dev = sqrt( sensors.mag.acc^2 + ...
    sensors.mag.lin^2 ) / sqrt(3) /100 * sensors.mag.FMR;  % standard deviation 
sensors.mag.variance = sensors.mag.STD_dev^2;              % variance 
sensors.mag.sensitivity = 100*1e3;                         % sensitivity [V/T]
sensors.mag.quant_V = 0.050;                                 % quantization interval [V]
sensors.mag.quant_T = sensors.mag.quant_V/sensors.mag.sensitivity;      % quantization interval [T]
sensors.mag.frequency = 5;                                    % update rate [Hz]
% magnetic sensor, non-orthogonality error
a = deg2rad(rand([3,1])*0.5); % orientation of the non othogonal reference system of the magnetometers [rad]
b = deg2rad(rand([3,1])*360);                           
sensors.mag.A_mag=[cos(a(1)) sin(a(1))*cos(b(1)) sin(a(1))*sin(b(1));... % rotation and distortion matrix
       sin(a(2))*cos(b(2)) cos(a(2)) sin(a(2))*sin(b(2));...
       sin(a(3))*cos(b(3)) sin(a(3))*sin(b(3)) cos(a(3))];    

%% actuators
magnetorquers.dipole = 350;                             % Maximum magnetic dipole [Am^2] da datsheet (attuatore su eoportal)

thruster.thrust = 0.01;                                 % Nominal thrust [N]   
thruster.direction = [ 0  0  0  0;   ...                % Direction in which the thrust is applied
                       0  0  0  0;   ...
                       1 -1  1 -1 ]; 
thruster.position = [ 0.3  0.3   -1.7    -1.7;   ...    % Thruster position wrt center of mass [m]
                     -1.05  1.05  1.1    -1.1;   ...
                      0.3  -0.3   0.195  -0.195  ];
thruster.firing_time = 0.05;                            % Minimum firing time [s]

%% Control
control.frequency = 20;                  %[Hz]
control.bdot_gain = 1e7*sat.I;
control.bdot_end = 0.01;
control.bdot_flag = false;
control.bdot_filter = 0.5;

control.pointing_toll = deg2rad(5);
control.pointing_Kp = -sat.I(3)/sat.I(1)*[1 1 1];
control.pointing_Kd = 0.5*sat.I .* [1 1 1]';

thruster.K_e = 47;
thruster.K_m = 40;
thruster.T_m = 20;
thruster.Uon = 19;
thruster.Uoff = 18.9;

