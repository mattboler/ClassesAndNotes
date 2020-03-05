clear all; close all; clc;

A = [-1 0; 0 1];

B = [0;1];
C = eye(2);
D = 0;

sys = ss(A, B, C, D);

x0 = [0;1];

t = [0:0.01:2];
c = -2.7183/1.7183;
u = c*ones(size(t));

lsim(sys, u, t, x0)