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

%% Choose algorithm and plots from command window

control.algorithm_vect = { 'De-tumbling', 'Pointing', 'De-tumbling + Pointing', 'No control' };

alg_index = input(['Select the algorithm you want to test typing its corresponding number: \n' ...
    ' 1: De-tumbling \n' ...
    ' 2: Slew and pointing \n' ...
    ' 3: De-tumbling, slew and pointing \n' ...
    ' 4: No control \n'],'s');
alg_index = str2double(alg_index);

if ~(alg_index == 1 || alg_index == 2 || alg_index == 3 || alg_index == 4)
    error('Choose a number between 1 and 4')
end

create_plots = input('Do you want to generate the plots? ("yes" or "no"): ','s');

if strcmp(create_plots,'yes')
    save_plots = input('Do you want to save the plots in .eps format? ("yes" or "no"): ','s');
end

configs;

%% Initial conditions after release from launcher

% initial angular velocity in body frame [rad/s]
settings.w0 = [3e-2,1e-2,3e-2]'; 

load_system('Model')
set_param('Model/Dynamics/Integrator2','InitialCondition','settings.w0') % set intial w in the simulink model

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

%% Run simulation

switch alg_index

    case 1 % De-tumbling

        control.algorithm = control.algorithm_vect{1};
        magnetorquers.dipole = 120;               % Maximum magnetic dipole [Am^2]
        phase1 = sim('Model.slx');

    case 2 % Pointing

        control.algorithm = control.algorithm_vect{2};
        settings.w0 = [1e-2 1e-2 1e-2]';
        magnetorquers.dipole = 350;      % Maximum magnetic dipole [Am^2]
        phase2 = sim('Model.slx');

    case 3 % De-tumbling + pointing

        %De-tumbling
        control.algorithm = control.algorithm_vect{1};
        magnetorquers.dipole = 120;               % Maximum magnetic dipole [Am^2]
        phase1 = sim('Model.slx');

        %Pointing
        % initial conditions ( final conditions of phase 1 )
        environment.theta_G0 = phase1.theta_G(end);
        r0 = phase1.r_N(end,:)';
        v0 = phase1.v(end,:)';
        settings.w0 = deg2rad(phase1.w_BN(end,:))';
        settings.E0 = deg2rad(phase1.E_312(end,:)');

        control.algorithm = control.algorithm_vect{2};
        magnetorquers.dipole = 350;      % Maximum magnetic dipole [Am^2]
        phase2 = sim('Model.slx');

    case 4 % No control
        
        Simtime = settings.Time;
        control.algorithm = control.algorithm_vect{4};
        magnetorquers.dipole = 0;               % Maximum magnetic dipole [Am^2]
        phase0 = sim('Model.slx');

end
%% plots
if strcmp(create_plots,'yes')
    plots;
end
