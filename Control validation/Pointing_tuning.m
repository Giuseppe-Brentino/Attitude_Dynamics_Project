
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
[E1,E2,E3] = dcm2angle(A_BN0,'ZXY');
A0 = angle2dcm(E1-0.5,E2+0.8,E3+0.3,'ZXY');

K_w = 1e-2;
K_cz = [1 1 1]; %sat.I/sat.I(3)*1e-4;
K_cx = [1 1 1]; %sat.I/sat.I(3)*1e-4;

thruster.K_e = 47;
thruster.K_m = 40;
thruster.T_m = 20;
thruster.Uon = 19;
thruster.Uoff = 18.9;

% sys = tf([thruster.K_m],[thruster.T_m 1]);
% sys = c2d(sys,0.05);
% num = sys.Numerator{1}
% den = sys.Denominator{1}

% %% GA
% ub = [100, 500]; % [K_m, T_m]
% lb = [0 0];
% 
% options = optimoptions('ga','PlotFcn','gaplotbestf','UseParallel', true, 'UseVectorized', false, 'PopulationSize', 100);
% X = ga(@(X)fitnessfcn(X,settings,sat,environment,thruster,n,r0,v0,A0,A_BN0),2,[],[],[],[],lb,ub,[],[],options);
% save Pointing_gains X