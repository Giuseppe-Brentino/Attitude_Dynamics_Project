clearvars; close all; clc;

load("WMM.mat")

N = 13; % ordine (N max is 12 because of K limit in WMM.mat)



R = 6371.2;
r = R+0;
lat = deg2rad(-80);
phi = deg2rad(240); %lambda longitudine

%lat = 30; % [deg] to test the code, it's an input in the simulink function
theta = (pi/2-lat); % theta is the colatitude [rad]
delta = pi/2-theta;
alpha_G = 0; %da cambiare su Simulink
alpha = deg2rad(phi + alpha_G);

%% Spherical Geocentric Coordinates
%test to validate the code

lambda = phi;
%phi = lat;
A = 6378137; %[m]
f = 1/298.257223563;
e = sqrt(f*(2-f));
R_c = A/sqrt(1-e^2*sin(lat)^2);
p = R_c*cos(lat);
z = sin(lat)*R_c*(1-e^2);
rr = sqrt(p^2+z^2);
theta = pi/2-asin(z/rr);



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

%% Magnetic Field

B_r = 0;
B_th = 0;
B_phi = 0;

%indice i scorre n da 1 a N
%indice j scorre m da 0 a n

for m = 0:N
    for n = 1:N
        if m<=n
            B_r = B_r + (R/r)^(n+2)*(n+1)*...
                ((g(n,m+1)*cos(m*phi) + h(n,m+1)*sin(m*phi))*P(n, m+1));
            B_th = B_th + (R/r)^(n+2)*...
                ((g(n,m+1)*cos(m*phi) + h(n,m+1)*sin(m*phi))*dP(n, m+1));
            B_phi = B_phi + (R/r)^(n+2)*...
                (m*(-g(n,m+1)*sin(m*phi) + h(n,m+1)*cos(m*phi))* P(n, m+1));

        end

    end
end

B_th = -B_th;
B_phi = -B_phi/sin(theta);

%Geocentric Coordinates
B_x = ( B_r*cos(delta) + B_th*sin(delta) ) * cos(alpha) - B_phi*sin(alpha);
B_y = ( B_r*cos(delta) + B_th*sin(delta) ) * sin(alpha) + B_phi*cos(alpha);
B_z = B_r*sin(delta) + B_th*cos(delta);

B = norm([B_r B_phi B_th]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO:
% Implementare K in modo da raggiungere ordini superiori a N=12
% Implementare il codice in modo più carino includendo P_01 e dP_01 nell'if (FACOLTATIVO)
% Aggiungere caso dei poli (singolarità)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%