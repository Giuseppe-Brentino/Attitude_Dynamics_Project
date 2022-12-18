clearvars; close all; clc;

addpath('..\functions\');
addpath('..\data\');
addpath('..\subsystems\');

configs;

%% de-tumbling
% linearizzazione intorno a w = [0, 0, 0]

dt = 0.1;                       % frequenza del controllore
w = 0;
I = sat.I;
settings.w0 = [0.0524 0.0524 0.0524]';
% initial attitude
h0=cross(r0,v0);
A_BN0 = [cross(-r0',h0')/norm(cross(r0,h0)); h0'/norm(h0); r0'/norm(r0)];
[E1, E2, E3] = dcm2angle(A_BN0, "ZXY");
settings.E0 = [E1 E2 E3]';

%% B_dot
k_b = 1e9;


