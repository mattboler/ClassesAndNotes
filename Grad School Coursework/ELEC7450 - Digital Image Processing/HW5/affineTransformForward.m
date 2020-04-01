function [outImage] = affineTransformForward(t_matrix, t_vector, t_image)
%AFFINETRANSFORM Perform an affine transform of an image using forward
%mapping

[h, w, d] = size(t_image);

% 1: Build arrays of coordinate vectors for easy transforming
inputVectors = getCoordinates(t_image);
%inputCoords = sub2ind(size(t_image), inputVectors(2,:), inputVectors(1,:));

% 1b: Shift vectors s.t. the center is the origin; this will force a
% rotation about the center of the input image.


du = floor(w/2);
dv = floor(h/2);
inputVectors_centered = inputVectors - [du; dv];

% 2: Apply transformation to vectors
outputVectors_centered = t_matrix * inputVectors_centered + t_vector;

% 3: Generate new image and use mapped vectors to populate image

% Get boundaries of output image (round towards making image larger to
% avoid cutting out information)
minOutU = floor(min(outputVectors_centered(1,:)));
minOutV = floor(min(outputVectors_centered(2,:)));
maxOutU = ceil(max(outputVectors_centered(1,:)));
maxOutV = ceil(max(outputVectors_centered(2,:)));

heightOut = maxOutV - minOutV + 1;
widthOut = maxOutU - minOutU + 1;

outImage = zeros(heightOut, widthOut, d);

% unshift output so top left is at (1,1)
padOutU = -1 * minOutU + 1; % If min val is 0: 1, if min val is 2: -1, if min val is -2: 3
padOutV = -1 * minOutV + 1;

outputVectors = outputVectors_centered + [padOutU; padOutV];
outputVectors = floor(outputVectors);
%outputCoords = sub2ind(size(outImage), outputCoords(2,:), outputCoords(1,:));

%outImage(outputCoords) = t_image(inputCoords);
inputVectors(:, 1:10);
outputVectors(:, 1:10);
for i = 1:size(inputVectors, 2)
    u_in = inputVectors(1,i);
    v_in = inputVectors(2,i);
    
    u_out = outputVectors(1,i);
    v_out = outputVectors(2,i);
    outImage(v_out, u_out, :) = t_image(v_in, u_in, :);
end

outImage = uint8(outImage)

end