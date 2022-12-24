%{
    Main script for the stochastic analysis of the de-tumbling time
%}

clearvars; close all; clc;

addpath('./subsystems/');
addpath('./data/');
addpath('./functions/');
addpath('./Control validation/');

configs;

Nsim = 400;         % number of simulations

%% Initial conditions after release from launcher

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

%% Create stochastic parameters

mu_w = 0;                                            % angular velocity's mean value
sigma_wx = (deg2rad(2))/3;                           % wx's standard deviation
sigma_wy = (deg2rad(3))/3;                           % wy's standard deviation
sigma_wz = (deg2rad(3))/3;                           % wz's standard deviation

% generate Nsim normally distributed initial angular velocities
stoch.wx = normrnd(mu_w,sigma_wx,1,Nsim);
stoch.wy = normrnd(mu_w,sigma_wy,1,Nsim);
stoch.wz = normrnd(mu_w,sigma_wz,1,Nsim);

control.algorithm = 'De-tumbling';
magnetorquers.dipole = 120;               % Maximum magnetic dipole [Am^2]
stoch.w0 = [stoch.wx;stoch.wy;stoch.wz];


%% define simulation variables
load_system('Model')
simIn(1:Nsim) = Simulink.SimulationInput('Model');
for i = 1:Nsim
    simIn(i) = simIn(i).setBlockParameter('Model/Dynamics/Integrator2',...
        'InitialCondition' , 'stoch.w0(:,i)');
    simIn(i) = simIn(i).setVariable('settings',settings);
    simIn(i) = simIn(i).setVariable('sat',sat);
    simIn(i) = simIn(i).setVariable('stoch',stoch);
    simIn(i) = simIn(i).setVariable('sensors',sensors);
    simIn(i) = simIn(i).setVariable('environment',environment);
    simIn(i) = simIn(i).setVariable('magnetorquers',magnetorquers);
    simIn(i) = simIn(i).setVariable('thruster',thruster);
    simIn(i) = simIn(i).setVariable('control',control);
    simIn(i) = simIn(i).setVariable('A_BN0',A_BN0);
    simIn(i) = simIn(i).setVariable('alfa1',alfa1);
    simIn(i) = simIn(i).setVariable('alfa2',alfa2);
    simIn(i) = simIn(i).setVariable('catalogue',catalogue);
    simIn(i) = simIn(i).setVariable('g',g);
    simIn(i) = simIn(i).setVariable('h',h);
    simIn(i) = simIn(i).setVariable('K',K);
    simIn(i) = simIn(i).setVariable('r0',r0);
    simIn(i) = simIn(i).setVariable('Simtime',Simtime);
    simIn(i) = simIn(i).setVariable('v0',v0);
    simIn(i) = simIn(i).setVariable('n',n);
    simIn(i) = simIn(i).setVariable('i',i);
end

%% Run simulation
simulation = parsim(simIn);
save simdata simulation
%% plots

stoch.wx = rad2deg(stoch.wx);
stoch.wy = rad2deg(stoch.wy);
stoch.wz = rad2deg(stoch.wz);

end_time = zeros(1,Nsim);
stoch.w = zeros(1,Nsim);
for i = 1:Nsim
    end_time(i) = simulation(1,i).tout(end)/3600;
    stoch.w(i) = norm([stoch.wx(i),stoch.wy(i),stoch.wz(i)]);
end

N_bar = min(Nsim,100); %number of bars in the histogram plot 

% time histogram
stoch_time_histo = figure('Units','normalized','OuterPosition',[0 0 1 1]);
histogram(end_time,N_bar)
xlabel('Time required for de-tumbling [h]')
ylabel('Number of occurrencies')

% w-t relationship

stoch_time_w = figure('Units','normalized','OuterPosition',[0 0 1 1]);

subplot(1,4,1)
fitting = fit( stoch.wx',end_time','poly2' );
plot(fitting,stoch.wx,end_time,'.')
xlabel('\omega_x [deg/s]')
ylabel('Time required for de-tumbling [h]')

subplot(1,4,2)
fitting = fit(stoch.wy',end_time','poly2' );
plot(fitting,stoch.wy,end_time,'.')
xlabel('\omega_y [deg/s]')
ylabel('Time required for de-tumbling [h]')

subplot(1,4,3)
fitting = fit( stoch.wz',end_time','poly2' );
plot(fitting,stoch.wz,end_time,'.')
xlabel('\omega_z [deg/s]')
ylabel('Time required for de-tumbling [h]')

subplot(1,4,4)
fitting = fit( stoch.w',end_time','poly2');
plot(fitting,stoch.w,end_time,'.')
xlabel('norm(\omega) [deg/s]')
ylabel('Time required for de-tumbling [h]')

% mean time and std
stoch_mean_std = figure('Units','normalized','OuterPosition',[0 0 1 1]);
mu_t = zeros(1,Nsim);
sigma_t = zeros(1,Nsim);
for i = 1:Nsim
    mu_t(i) = mean(end_time(1:i));
    sigma_t(i) = std(end_time(1:i));
end

subplot(1,2,1)
plot(1:Nsim,mu_t)
xlabel('Number of iterations')
ylabel("Stop time's mean value [h]")
subplot(1,2,2)
plot(1:Nsim,sigma_t)
xlabel('Number of iterations')
ylabel("Stop time's standard deviation [h]")

exportgraphics(stoch_time_w,'.\report_images\stoch_time_w.pdf','ContentType','vector')
exportgraphics(stoch_time_histo,'.\report_images\stoch_time_histo.pdf','ContentType','vector')
exportgraphics(stoch_mean_std,'.\report_images\stoch_mean_std.pdf','ContentType','vector')