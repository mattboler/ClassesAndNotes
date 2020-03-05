clc; clear all; close all;

v_drone = 100; % km/hr
R_0 = 6378137;
circumference = 2 * pi * R_0;

time = circumference / (v_drone * 10^3);