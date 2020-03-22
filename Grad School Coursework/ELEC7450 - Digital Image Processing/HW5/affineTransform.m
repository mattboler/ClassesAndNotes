function [outImage] = affineTransform(t_matrix, t_vector, t_image)
%AFFINETRANSFORM Summary of this function goes here
%   Detailed explanation goes here



% 1: Build arrays of coordinates for easy transforming
input_coords = getCoordinates(t_image);


end

function coords = getCoordinates(t_input)
% Get coordinates in [x1 x2 x3 ...; y1 y2 y3 ...] form
    [h, w] = size(t_input);
    [U, V] = meshgrid(1:w, 1:h);
    U_row = reshape(U, 1, numel(U));
    V_row = reshape(V, 1, numel(V));
    coords = [U_row; V_row];
end

function [h_out, w_out] = getOutputSize(t_matrix, t_vector, t_image)
    [h, w] = size(t_image);
    extreme_coords = [0 0 w w; 0 h 0 h];
    out_coords = t_matrix * extreme_coords + t_vector;
    h_out = max(out_coords(2,:)) - min(out_coords(2,:));
    w_out = max(out_coords(1,:)) - min(out_coords(1,:));
end