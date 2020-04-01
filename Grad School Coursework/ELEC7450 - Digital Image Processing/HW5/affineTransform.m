function [outImage] = affineTransform(t_matrix, t_vector, t_image)
%AFFINETRANSFORM Summary of this function goes here
%   Detailed explanation goes here

% Treating center of images as origin -> Need to shift to center before rotating!
% We transform the extreme points of the input to find the boundaries of
% the output
inputCoords = getCoordinates(t_image);
[outputCoords, outImage] = findOutputCoords(t_matrix, t_vector, t_image)

% Use homogeneous transform to make inversion simple
T = [t_matrix, t_vector; 0 0 1];
T_inv = inv(T);
A_back = T_inv(1:2, 1:2);
b_back = T_inv(1:2,3);

backwarpedOutputCoords = A_back * (outputCoords - (size(outImage)/2)')  + b_back + (size(outImage)/2)'
end

function coords = getCoordinates(t_input)
% Get coordinates in [x1 x2 x3 ...; y1 y2 y3 ...] form
    [h, w] = size(t_input);
    [U, V] = meshgrid(1:w, 1:h);
    U_row = reshape(U, 1, numel(U));
    V_row = reshape(V, 1, numel(V));
    coords = [U_row; V_row];
end

function [coords, outImage] = findOutputCoords(t_matrix, t_vector, t_image)
[h, w] = size(t_image);
extremePoints = [0 0 w w; 0 h 0 h];
% Pre-shift to move origin to center
extremePoints = extremePoints - [w/2; h/2];

% Apply transform
extremePointsTransformed = t_matrix * extremePoints + t_vector

% Shift back
extremePointsTransformed = extremePointsTransformed + [w/2; h/2];

height = ceil( max(extremePointsTransformed(2,:)) - min(extremePointsTransformed(2,:)) );
width = ceil( max(extremePointsTransformed(1,:)) - min(extremePointsTransformed(1,:)) );

outImage = zeros(height, width);
coords = getCoordinates(outImage);

end