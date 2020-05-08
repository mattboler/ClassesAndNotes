clc; clear all; close all;

%% Part 1: LQR + KF

dt = 0.001;

% State Vector : [x, v]

m = 10; % kg

v_sigma = 0.1; % m
w_sigma = 1.0; % N

A = [0 1; 0 0];
B = [0; 1/m];
B_w = [0; 1/m];
C = [1 0];
D = 0;

x_0 = [10; 0];

