function [projected_points] = project_points(points, K, R, t)
%PROJECT_POINTS Summary of this function goes here
%   points : 3xn nonhomogeneous array
%   K : 3x3 camera intrinsic matrix
%   R : 3x3 rotation matrix
%   t : 3x1 translation vector

pose = [R, t];
projection_matrix = K * pose;

[dim_points, num_points] = size(points);

w_component = ones(1, num_points);
temp_points = [points; w_component];

unscaled_projected_points = projection_matrix * temp_points;

projected_points = zeros(dim_points - 1, num_points);

projected_points(1,:) = unscaled_projected_points(1,:) ./ unscaled_projected_points(3,:);
projected_points(2,:) = unscaled_projected_points(2,:) ./ unscaled_projected_points(3,:);


end

