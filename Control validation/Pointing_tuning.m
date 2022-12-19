
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
K_cz = sat.I/sat.I(3)*1e-4;%[1 1 1];
K_cx = sat.I/sat.I(3)*1e-4;%[1 1 1];

settings.E0 = [E1 E2 E3]';                % initial euler angles [rad]

%% GA
ub = 1000*ones(7,1);
lb = zeros(7,1);

options = optimoptions('ga','PlotFcn','gaplotbestf','UseParallel', true, 'UseVectorized', false);
X = ga(@(X)fitnessfcn(X,settings,sat,environment,n,r0,v0,A0,A_BN0),7,[],[],[],[],lb,ub,[],[],options);
save Pointing_gains X
