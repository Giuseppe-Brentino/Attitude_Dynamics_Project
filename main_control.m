clearvars; close all; clc;

addpath('./subsystems/');
addpath('./data/');
addpath('./functions/');


main;

%% LQR de-tumbling
% linearizzazione intorno a w = [0, 0, 0]

dt = 0.1;                       % frequenza del controllore
A = eye(3);
B = dt*inv(diag(sat.I));

