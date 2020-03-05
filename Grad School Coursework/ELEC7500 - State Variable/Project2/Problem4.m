clc; clear all; close all;

A = [ 0 1; 2 1];

[V, D] = eig(A);

S = D;
M = V;
syms t;

phi = M*expm(S*t)*inv(M);
