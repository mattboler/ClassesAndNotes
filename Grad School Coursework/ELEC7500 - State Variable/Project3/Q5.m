clear all; close all; clc;

A = [1 0 0 0; 0 -2 0 0; 0 0 3 0; 0 0 0 -4];
B = [1;0;1;0];

syms k1 k3;
K = [k1 0 k3 0];

Acl = A - B*K;

eigs = eig(Acl);
eqn1 = eigs(3) == -5;
eqn2 = eigs(4) == -5;



Y = solve([eqn1 eqn2], [k1 k3])

