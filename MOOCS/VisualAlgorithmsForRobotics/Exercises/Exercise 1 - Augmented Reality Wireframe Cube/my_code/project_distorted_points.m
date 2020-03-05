function projected_points = project_distorted_points(points, K, R, t, D)
%PROJECT_DISTORTED_POINTS Summary of this function goes here
%   Detailed explanation goes here
pose = [R, t];

[dim_points, num_points] = size(points);

w_component = ones(1, num_points);
temp_points = [points; w_component];

camera_points = pose * temp_points;

% Get normalized points (not including centering)
camera_points_normalized = [camera_points(1,:) ./ camera_points(3,:);
    camera_points(2,:) ./ camera_points(3,:)];

% Apply distortion
k1 = D(1);
k2 = D(2);

r2 = camera_points_normalized(1,:).^2 + camera_points_normalized(2,:).^2;
r4 = r2.^2;
adder = ones(size(r2));

coefficient_vector = adder + k1*r2 + k2*r4;
coefficient_matrix = [coefficient_vector; coefficient_vector];

prime_coords = [camera_points_normalized .* coefficient_matrix; adder];

lambda_coords = K * prime_coords;

projected_points = [lambda_coords(1,:) ./ lambda_coords(3,:);
    lambda_coords(2,:) ./ lambda_coords(3,:)];

end

