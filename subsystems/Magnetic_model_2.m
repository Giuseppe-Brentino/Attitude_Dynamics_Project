clearvars; close all; clc;

N = 13; % ordine 

n = 1:13;
m = 0:13;

theta = deg2rad(30); % to test the code, it's an input in the simulink function


%% P-dP
% uguale a fcn su matlab ma non a paper

P = zeros(N,N+1);
dP = zeros(N,N+1);

P_00 = 1;
dP_00 = 0;

P(1,1) = cos(theta);
dP(1,1) = -sin(theta)*P_00;

P(1,2) = sin(theta);  
dP(1,2) = cos(theta)*P_00;

for i = 2:13
    for j = 0:i
        if j == i
        P(i,j+1) = sin(theta)*P(i-1,j);
        dP(i,j+1) = sin(theta)*dP(i-1,j) + cos(theta)*P(i-1,j);
        elseif i < 3
            P(i,j+1) = cos(theta)*P(i-1,j+1);
            dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1);
        else
            P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P(i-2,j+1);
            dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1) ...
                - K(i,j+1)*dP(i-2,j+1);
        end
    end
end

