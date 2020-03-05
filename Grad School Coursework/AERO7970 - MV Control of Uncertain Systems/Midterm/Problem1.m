clc; clear all; close all;

A = [-0.2, -100, 0, 0, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 0;
    0, 0, -.1, -10, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 0, 0;
    0, 0, 0, 0, 0, -5, -6;
    0, 0, 0, 0, 0, 1, 0];

B = [1, 0;
    0, 0;
    1, 0;
    0, 0;
    0, 1;
    0, 1;
    0, 0];

C = [10, 10, 0, 0, 1, 0, 0;
    0, 0, 1, 2, 0, 5, 5];

D = 0;

sys = ss(A, B, C, D);
G = pck(sys.A, sys.B, sys.C, sys.D);

CT = rank(ctrb(A, B))
OB = rank(obsv(A, C))

trans_zeros = tzero(sys)

eigs = eig(A)

h2 = norm(sys)

hinf = norm(sys, Inf)

sigma(G),grid, legend

%{
w = logspace(-1,1,200);
Gf = frsp(G, w);
[u,s,v] = vsvd(Gf);
vplot('liv,lm', s), grid
%}