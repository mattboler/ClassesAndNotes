clc; clear all; close all;

% Portland coordinates
port_lat = 45; % deg
port_lon = -120; % deg

% Tashkent coordinates
tash_lat = 45; % deg
tash_lon = 60; % deg

% Assumptions:
% 10km altitude
% 800km/hr speed
% Spherical earth
R_0 = 6378137;
alt = 10 * 10^3;
lat = 45;
speed = 800;

% Try going around?

R_eff = (R_0 + alt) * cosd(lat);

circumference = pi * (2*R_eff);

dist = (tash_lon - port_lon)/360 * circumference;
time = dist / (speed * 10^3);

% But wait, over the top of the earth is only 90 degrees rather than 180!
% Assuming nonrotating earth
R_eff2 = R_0 + alt;
dist2 = (port_lat + tash_lat)/360 * 2 *R_eff2 * pi;
time2 = dist2 / (speed * 10^3);


