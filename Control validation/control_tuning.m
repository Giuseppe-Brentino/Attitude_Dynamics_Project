
clearvars; close all; clc;

addpath('../subsystems/');
addpath('../data/');
addpath('../functions/');


configs;

%% Nadir pointing
w0 = [1.5e-3 1.5e-3 1.5e-3]';                % initial angular velocity in body frame [rad/s]

% initial attitude
A_BN0 = [v0'/norm(v0); (cross(r0, v0))'/norm(cross(r0, v0)); r0'/norm(r0)];

n = sqrt(settings.mu/settings.a^3);
A = [                    0            0   ( sat.I(2)-sat.I(3) )/sat.I(1)*n;...
                         0            0                       0;...
     ( sat.I(1)-sat.I(2) )/sat.I(3)*n 0                       0 ];

C = eye(3);
B = diag([1/sat.I(1) 1/sat.I(2) 1/sat.I(3)]);
D = zeros(3);

% lb = - 1*ones(6,1);
% ub =   1*ones(6,1);
% options = optimoptions("ga","PlotFcn","gaplotbestf");
% X = ga(@(X)fitnessfcn(X,n,A,B,C,D,settings,r0,v0,A_BN0,w0),6,[],[],[],[],lb,ub,[],options);

% % load("Nadir.mat")
% % Kp = 7*X(1:3);
% % Kd = 2*X(4:6);

Kp = -[-0.1 -0.1 0]';
Kd = -[10 10 0];
