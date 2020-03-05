function R = extract_rotation(pts_1,pts_2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
focalLength = 400; % arbitrary value, only adjusts z value in orthogonal case

% Determine Camera Paramters
[rows1, columns1] = size(I_t1);
IntrinsicMatrix = [focalLength,0,0;0,focalLength,0;columns1/2,rows1/2,1];
radialDistortion = [0,0,0]; 
cameraParams = cameraParameters('IntrinsicMatrix',IntrinsicMatrix,'RadialDistortion',radialDistortion);

[tform, inlier1, inlier2] = estimateGeometricTransform(pts_1, pts_2, 'affine');
[relOr, relLoc] = relativeCameraPose(tform, cameraParams, inlier1, inlier2);
R = relOr;

end

