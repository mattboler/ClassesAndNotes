clc; clear all; close all;

path = '/home/matt/Projects/VisualAlgorithms/Exercises/Exercise 1 - Augmented Reality Wireframe Cube';

K = dlmread([path '/data/K.txt']);
D = dlmread([path '/data/D.txt']);
poses = dlmread([path '/data/poses.txt']);

%% Build world coords of checkerboard
% (0,0,0) is the upper left corner
x_coords = 0:8;
y_coords = 0:5;
scale = 0.04; % cm

[x, y] = meshgrid(x_coords * scale, y_coords * scale);

[m, n] = size(x);

n_points = m * n;

world_points = zeros(3, n_points);

for i = 1:n_points
   % Using linear coordinates for simplicity
   world_points(1,i) = x(i);
   world_points(2,i) = y(i);
end


%% Undistorted image
img = imread('data/images_undistorted/img_0001.jpg');
pose = poses(1,:);
rvec = pose(1:3);
t = pose(4:6)';

r_mat = convert_rvec_to_matrix(rvec);
image_points = project_points(world_points, K, r_mat, t);

u = image_points(1,:);
v = image_points(2,:);


%% Distorted image
img = imread('data/images/img_0001.jpg');
image_points = project_distorted_points(world_points, K, r_mat, t, D);
u = image_points(1,:);
v = image_points(2,:);

figure(2);
imshow(img);
hold on;
scatter(u, v, 'b*');
hold off;



