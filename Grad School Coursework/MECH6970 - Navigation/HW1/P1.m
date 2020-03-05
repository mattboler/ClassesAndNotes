clc; clear all; close all;

p_ecef = [-4.2198; -26.428; 9.4295]*10^6; % m
v_ecef = [2941.4; 483.43; 2667.2]; % m/s

%% A

% What are the satellite's geodetic position coordinates?

R_0 = 6378137;
E = 0.0818191908426215;

x = p_ecef(1);
y = p_ecef(2);
z = p_ecef(3);

lon = atan2d(y, x);

prev_lat = 0;
lat = 10;

count = 0;

while abs(lat - prev_lat) > 0.00001
    count = count + 1;
    prev_lat = lat;
    N = (R_0)/sqrt( 1 - (E*sind(lat))^2 );
    H = -N + ( sqrt( x^2 + y^2 )/cosd(lat) );
    lat = atan2d( (z / sqrt( x^2 + y^2 )) , ( 1 - (E^2)*(N / (N+H)) ));
end

%% B 

% No computations needed

%% C

% What is the speed in km/s?

speed_m = norm(v_ecef);
speed_km = speed_m * 10^(-3);

%% D

% Velocity in a NED system below the satellite?

R1 = [cosd(lon), sind(lon), 0; -sind(lon), cosd(lon), 0; 0, 0, 1];
alpha = lat + 90;
R2 = [cosd(alpha), 0, sind(alpha); 0, 1, 0; -sind(alpha), 0, cosd(alpha)];

C_e_ned = R2 * R1;
v_ned = C_e_ned * v_ecef;

%% E

% % Az/El for an observer at [0 8 28; 78 29 17; 9228]?
o_lat = [-0 -8 -28];
o_lon = [-78 -29 -17];
o_h = 9228/3.281;

o_lat = o_lat(1) + o_lat(2)/60 + o_lat(3)/3600;
o_lon = o_lon(1) + o_lon(2)/60 + o_lon(3)/2600;
% 
% obsv_lla = [o_lat; o_lon; o_h];
% sat_lla = [lat; lon; H];
% 
R1_obsv = [cosd(o_lon), sind(o_lon), 0; -sind(o_lon), cosd(o_lon), 0; 0, 0, 1];
alpha_obsv = o_lat + 90;
R2_obsv = [cosd(alpha_obsv), 0, sind(alpha_obsv); 0, 1, 0; -sind(alpha_obsv), 0, cosd(alpha_obsv)];

C_obsv = R2_obsv * R1_obsv;
% obsv_ned = [0;0;-9228 * .3048];
% sat_ned = C_obsv * p_ecef;
% 
% del_ned = sat_ned - obsv_ned;
% 
% az = atan2d(del_ned(2), del_ned(1));
% rise = -del_ned(3);
% run = sqrt(del_ned(1)^2 + del_ned(2)^2);
% el = atan2d(rise, run);

Re = R_0 / (sqrt(1 - (E^2)*sind(o_lat)^2));
p_ECEF_Quito = [(Re+o_h)*cosd(o_lat)*cosd(o_lon), (Re+o_h)*cosd(o_lat)*sind(o_lon), ...
    (Re*(1-E^2)+o_h)*sind(o_lat)]';

relative_p_ECEF = p_ecef - p_ECEF_Quito;
p_NED = C_obsv * relative_p_ECEF;
range = norm(p_NED);
azimuth = rad2deg(atan2(p_NED(2), p_NED(1)));
elevation = rad2deg(asin(-p_NED(3)/range));

