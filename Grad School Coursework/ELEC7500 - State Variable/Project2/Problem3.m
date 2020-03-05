clc; clear all; close all;

A = [2 -1; 3 4];
syms s;
I = eye(2);

char_mat = s*I - A;
char_poly = det(char_mat);

disp("A^5")
-924*eye(2) + 229*A
A^5

disp("A^-1")
(6/11)*eye(2) - (1/11)*A
inv(A)

disp("ln(A)");

