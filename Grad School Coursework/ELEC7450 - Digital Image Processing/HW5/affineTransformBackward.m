function [outImage] = affineTransformBackward(t_matrix, t_vector, t_image)
%AFFINETRANSFORMBACKWARD Summary of this function goes here
%   Detailed explanation goes here

[h, w, d] = size(t_image);
du_in = floor(w/2);
dv_in = floor(h/2);

T = [t_matrix, t_vector; 0 0 1];
T_inv = inv(T);
A_back = T_inv(1:2, 1:2);
b_back = T_inv(1:2, 3);

% Part 1: Generate empy output image
% Transform corners of input image to find bounds of output image

du = floor(w/2);
dv = floor(h/2);
inputCorners = [0 0 w w; 0 h 0 h] - [du; dv];
outputCorners = t_matrix * inputCorners + t_vector;

minOutU = floor(min(outputCorners(1,:)));
minOutV = floor(min(outputCorners(2,:)));
maxOutU = ceil(max(outputCorners(1,:)));
maxOutV = ceil(max(outputCorners(2,:)));

% This is the size of the area covered by the transformed input image
h_out = maxOutV - minOutV + 1;
w_out = maxOutU - minOutU + 1;
outImage = zeros(h_out + ceil(t_vector(2)), w_out + ceil(t_vector(1)), d); 

% Part 2: Pad input image with reflected copies of original
% Backwards warp corners of output to find where/if they exceed input image

zero_image = zeros(size(t_image));

top_bot = zero_image;
left_right = zero_image;
corners = zero_image;
padded_image = [corners, top_bot, corners; left_right, t_image, left_right; corners, top_bot, corners];

% Part 3: Backwards map!

outCoords = getCoordinates(outImage);
du_out = floor(w_out/2);
dv_out = floor(h_out/2);
outCoords_centered = outCoords - [du_out; dv_out];

for i = 1:size(outCoords_centered, 2)
   u = outCoords_centered(1,i);
   v = outCoords_centered(2,i);
   
   p = [u; v];
   p_inv = A_back * p + b_back;
   % Coords in uncentered input image
   u_in = round(p_inv(1) + du_in + w);
   v_in = round(p_inv(2) + dv_in + h);
   
   outImage(v+dv_out, u+du_out, :) = padded_image(v_in, u_in, :);
end

outImage = uint8(outImage);

end

function coords = getCoordinates(t_input)
% This function generates a [u; v] coordinate vector for each pixel in an
% image
    [h, w, ~] = size(t_input);
    [U, V] = meshgrid(1:w, 1:h);
    U_row = reshape(U, 1, numel(U));
    V_row = reshape(V, 1, numel(V));
    coords = [U_row; V_row];
end
