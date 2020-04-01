function coords = getCoordinates(t_input)
% This function generates a [u; v] coordinate vector for each pixel in an
% image
    [h, w, d] = size(t_input);
    [U, V] = meshgrid(1:w, 1:h);
    U_row = reshape(U, 1, numel(U));
    V_row = reshape(V, 1, numel(V));
    coords = [U_row; V_row];
end
