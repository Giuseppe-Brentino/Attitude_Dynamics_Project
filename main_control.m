clearvars; close all; clc;

addpath('./subsystems/');
addpath('./data/');
addpath('./functions/');


main;

%% LQR de-tumbling
% linearizzazione intorno a w = [0, 0, 0]

dt = 0.1;                       % frequenza del controllore
w = 0;
I = sat.I;


%% test ga
% parpool
% A =[0 (I(2)-I(3))/I(1)*w (I(2)-I(3))/I(1)*w;...
%     (I(3)-I(1))/I(2)*w 0 (I(3)-I(1))/I(2)*w;...
%     (I(1)-I(2))/I(3)*w (I(1)-I(2))/I(3)*w 0] ;
% final_w = 0.002*ones(3,1); % 0.1 deg/s
% high_torque = 0.03; %Nm
% Ts = 2800;
% B = inv(diag(I));
% lb = [0 1 1];
% ub = [1 100 100];
% options = optimoptions("ga","PlotFcn","gaplotbestf");
% X = ga(@(X)fitnessfcn(X,A,B,final_w,high_torque,Ts,sat,settings),3,[],[],[],[],lb,ub,[],options);

%% LQR aaaa
%  load("LQR_Data4.mat")
% % A =[0 (I(2)-I(3))/I(1)*w (I(2)-I(3))/I(1)*w;...
% %     (I(3)-I(1))/I(2)*w 0 (I(3)-I(1))/I(2)*w;...
% %     (I(1)-I(2))/I(3)*w (I(1)-I(2))/I(3)*w 0] + X(1)*eye(3);
% % B = inv(diag(sat.I));
% % R = diag(X(3)*ones(3,1));
% % Q = diag(X(2)*ones(3,1));
% % 
% % [k,s,clp] = lqr(A,B,Q,R);
% % 
 save LQR_Data5 X

%% LQR AAAAAAAAAAAAAAAAAAAAA

% A =[0 (I(2)-I(3))/I(1)*w (I(2)-I(3))/I(1)*w;...
%     (I(3)-I(1))/I(2)*w 0 (I(3)-I(1))/I(2)*w;...
%     (I(1)-I(2))/I(3)*w (I(1)-I(2))/I(3)*w 0] +  X(1)*eye(3);%diag([0.2 0.4 0.2]);
% B = inv(diag(sat.I));
% R = diag(X(3)*ones(3,1));
% Q = diag(X(2)*ones(3,1));
% 
% [k,s,clp] = lqr(A,B,Q,R);
% 
% save k2.mat k;