
clearvars; close all; clc;

addpath('../subsystems/');
addpath('../data/');
addpath('../functions/');


configs;

%% Nadir pointing
settings.w0 = [1e-2 1e-2 1e-2]';                % initial angular velocity in body frame [rad/s]

% initial attitude
A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];
n = sqrt(settings.mu/settings.a^3);
[E1,E2,E3] = dcm2angle(A_BN0);
A0 = angle2dcm(E1+0.7,E2+0.4,E3-1);
% A0 = A_BN0;

K_w = 1e-2;
K_cz = [1 1 1]; %sat.I/sat.I(3)*1e-4;
K_cx = [1 1 1]; %sat.I/sat.I(3)*1e-4;

settings.E0 = [E1 E2 E3]';                % initial euler angles [rad]

%% Thruster control
thruster.K_e = 0.1;
thruster.K_m = 10;
thruster.T_m = 40;
thruster.Uon = 15;
thruster.Uoff = 5;


% %% GA
% ub = [100, 500]; % [K_m, T_m]
% lb = [0 0];
% 
% options = optimoptions('ga','PlotFcn','gaplotbestf','UseParallel', true, 'UseVectorized', false, 'PopulationSize', 100);
% X = ga(@(X)fitnessfcn(X,settings,sat,environment,thruster,n,r0,v0,A0,A_BN0),2,[],[],[],[],lb,ub,[],[],options);
% save Pointing_gains X