clc; clear all; close all;

%% Part 1:
A_cl = [0 1; -1 -1.4];
C = [1 0];
B_w = [0; 1];

dist_sigma = 2;
noise_sigma = 1;

dt = 0.1;
Q = ones(2) * dist_sigma * dist_sigma; Q(1,1) = 0; % equivalent to B_w * Q * B_w
R = ones(1) * noise_sigma * noise_sigma;

Q_d = B_w' * Q * B_w * dt;
A_cl_d = expm(A_cl * dt);

x(:,1) = [0; 0];

for i = 2 : 100/dt
    sensor_noise = randn * noise_sigma;
    position_noise = randn * dist_sigma;
    velocity_noise = randn * noise_sigma;
    
    corrupted_state = x(:, i-1) + 
    
    x(i) = A_cl_d * x(:, i-1) + B_w * [position_noise; velocity_noise];
    
    
end