clearvars; close all; clc;

load("WMM.mat")

N = 12; % ordine (N max is 12 because of K limit in WMM.mat)

n = 1:13;
m = 0:13;

lat = 30; % [deg] to test the code, it's an input in the simulink function
theta = deg2rad(90-lat); % theta is the colatitude [rad]

%% P-dP
% validated through the P and dP tables in "Jeremy Davis, Mathematical Modeling of Earthas Magnetic
% Field, 2004"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_{(n),(m+1)}: 
% given P(a,b) -> n=a, m=b-1
% n from 1 to N
% m from 0 to N
% N order of accuracy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(N,N+1);
dP = zeros(N,N+1);

P_00 = 1;
dP_00 = 0;

P_01 = 0; %?? ho messo questo valore perchè è quello che sta nelle tabelle del pdf ma ha senso?
dP_01 = 0;

P(1,1) = cos(theta); % n=1, m=0
dP(1,1) = -sin(theta)*P_00; 

P(1,2) = sin(theta);  % n=1, m=1
dP(1,2) = cos(theta)*P_00;

for i = 2:N % ciclo in n a partire da n=2 a n=N
    for j = 0:i % ciclo in m da m=0 a m=n
        if j == i
        P(i,j+1) = sin(theta)*P(i-1,j);
        dP(i,j+1) = sin(theta)*dP(i-1,j) + cos(theta)*P(i-1,j);
        elseif i < 3
            if j == 0 % n=2, m=0
            P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P_00;
            dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1) ...
                - K(i,j+1)*dP_00;
            elseif j == 1 % n=2, m=1
                P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P_01;
                dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1)...
                    - K(i,j+1)*dP_01;
            end
        else
            P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P(i-2,j+1);
            dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1) ...
                - K(i,j+1)*dP(i-2,j+1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO:
% Implementare K in modo da raggiungere ordini superiori a N=12
% Implementare il codice in modo più carino includendo P_01 e dP_01 nell'if (FACOLTATIVO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%