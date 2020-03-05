clc, clear all, close all;

%% P1:
% A: Calculate err between toolbox and class estimate for camera
% translation and rotation between Google-maps-ref.jpg and cropped image of
% the same.

num_features = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000];
[calculated_orientations, simulated_orientations, calculated_locations, simulated_locations] = generateData();

translation_err = [];
rotation_err = [];

for i = 1:length(num_features)
    translation_err(i) = norm(calculated_locations{i} - simulated_locations{i});
end

for i = 1:length(num_features)
   R1 = calculated_orientations{i};
   R2 = simulated_orientations{i};
   diff_R = R1 * R2';
   r = dcm2rod(diff_R);
   rotation_err(i) = norm(r);
end

% B: Plot a graph of these two errors as a function of number of features
% matched between the two images
figure(1);
plot(num_features, translation_err);
title("Translation error as a function of features matched");

figure(2);
plot(num_features, rotation_err);
title("Rotation error as a function of features matched");

%% P2: 
% Estimate the rotation matrix against Shelby-ref.jpg
% Can you describe what happened to the camera in each case?