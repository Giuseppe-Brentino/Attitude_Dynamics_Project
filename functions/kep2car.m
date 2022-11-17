function [r, v] = kep2car(settings)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Conversion from Keplerian elements to Cartesian coordinates. Angles in
%  radians
%  
% INPUT:
%   a   [1x1] Semi-major axis           [km]
%   e   [1x1] Eccentricity              [-]
%   i   [1x1] Inclination               [rad]
%   OM  [1x1] RAAN                      [rad]
%   om  [1x1] Pericentre anomaly        [rad]
%   th  [1x1] True anomaly              [rad]
%   mu  [1x1] Gravitational parameter   [km^3/s^2]
%
% OUTPUT:
%   r   [3x1] Position vector           [km]
%   v   [3x1] Velocity vector           [km/s]
%
% Contributors: 
%   Nicolò Galletta, Virginia di Biagio Missaglia
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = settings.mu;
a = settings.a;
e = settings.e;
i = settings.i;
OM = settings.OM;
om = settings.om;
th = settings.theta;


p = a*(1-e^2);          % semi-latus rectum [km]
r = p/(1+e*cos(th));    % radius [km]

% State vector in perifocal reference frame
r_PF = r*[cos(th), sin(th), 0]';
v_PF = sqrt(mu/p)*[-sin(th), e+cos(th), 0]';


%% Rotation of the state vector in the Earth Centered Equatorial Inertial Frame

%Starting frame: ECEI
% Rotation around axis k of an angle OM
R3_OM=[cos(OM) sin(OM) 0;
      -sin(OM) cos(OM) 0;
       0       0       1];

% Rotation around i'= N (node line) of an angle i
R1_i=[1  0      0;
      0  cos(i) sin(i);
      0 -sin(i) cos(i)];

% Rotation around k'' = h (angular momentum vector) of an angle om
R3_om=[ cos(om) sin(om) 0;
       -sin(om) cos(om) 0;
        0       0       1];

% Final frame: Perifocal frame
T_PF2ECI=[R3_om*R1_i*R3_OM]';     % Complete rotation matrix from PF to ECI

r = T_PF2ECI * r_PF;      %r_ECI
v = T_PF2ECI * v_PF;      %v_ECI

return




