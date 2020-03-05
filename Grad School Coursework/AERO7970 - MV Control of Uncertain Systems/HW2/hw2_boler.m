clc, clear all;

%{
AERO 7970 HW2
Matt Boler
6/19/2019
%}

%{
Problem 1

Given the following system:
%}
A_1 = [0,1,0,0,0;
       0,0,1,0,0;
       0,0,0,1,0;
       0,0,0,0,1;
       0,0,0,0,0];
[n_1, ~] = size(A_1);
B_1 = [0;
       0;
       0;
       0;
       1];

%{
1. Check the controllability of the system
2. Check the stability of the system
3. What is the response to initial conditions?
%}
       
C = ctrb(A_1, B_1);
disp("Controllability of Problem 1:");
fprintf("\tDim of A = %d, Rank of C = %d \n", n_1, rank(C));

% Note: Should I be doing more here? 
% Since it's controllable, I can just place poles wherever
lambda_1 = eig(A_1);
disp("Stability of Problem 1:");
fprintf("\tEigenvalues of A: %d %d %d %d %d \n", lambda_1);
fprintf("\tNumber of eigenvalues with positive real components: %d \n", sum(real(lambda_1) > 0));
fprintf("\tNumber of eigenvalues with negative real components: %d \n", sum(real(lambda_1) < 0));
fprintf("\tNumber of eigenvalues with 0 real component: %d \n", sum(real(lambda_1) == 0));

% TODO: Response to initial conditions
syms t;
phi_1 = expm(A_1*t);

%{
Problem 2

Given the following system:
%}
A_2 = [0,1,0,0,0;
       0,0,1,0,0;
       0,0,0,1,0;
       0,0,0,0,1;
       0,0,0,0,0];
[n_2, ~] = size(A_2);
C_2 = [1,0,0,0,0];

%{
1. Check the observability of the system
2. Check the stability of the system
3. What is the response to initial conditions?
%}

O = obsv(A_2, C_2);
disp("Observability of Problem 2:");
fprintf("\tDim of A = %d, Rank of O = %d \n", n_2, rank(O));

%{
Problem 3

Given the following system:
%}

syms w_0 r_0 m
A_3 = [0,       1,          0, 0,         0,      0;
       3*w_0^2, 0,          0, 2*w_0*r_0, 0,      0;
       0,       0,          0, 1,         0,      0;
       0,       -2*w_0/r_0, 0, 0,         0,      0;
       0,       0,          0, 0,         0,      1;
       0,       0,          0, 0,         -w_0^2, 0];

B_3 = [0,   0,         0;
       1/m, 0,         0;
       0,   0,         0;
       0,   1/(m*r_0), 0;
       0,   0,         0;
       0,   0,         1/(m*r_0)];
   
C_3 = [1, 0, 0, 0, 0, 0;
       0, 0, 1, 0, 0, 0;
       0, 0, 0, 0, 1, 0];
   
D_3 = 0;

[n_3, ~] = size(A_3);

%{
Determine the stability, controllability, and observability properties
%}

lambda_3 = eig(A_3);

disp("Stability of Problem 3:");
fprintf("\tEigenvalues of A: %s %s %s %s %s %s \n", char(lambda_3));
fprintf("\n");
beta_3 = [B_3, A_3*B_3, A_3^2 * B_3, A_3^3 * B_3, A_3^4 * B_3, A_3^5 * B_3];

disp("Controllability of Problem 3:");
fprintf("\tDim of A = %d, Rank of C = %d \n", n_3, rank(beta_3));

gamma_3 = [C_3;
           C_3*A_3;
           C_3*A_3^2;
           C_3*A_3^3;
           C_3*A_3^4;
           C_3*A_3^5];

disp("Observability of Problem 3:");
fprintf("\tDim of A = %d, Rank of O = %d \n", n_3, rank(gamma_3));   
       
       

