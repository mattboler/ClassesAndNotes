function [rotation_matrix] = convert_rvec_to_matrix(rvec)
%CONVERT_RVEC_TO_MATRIX Rodrigues formula to convert axis-angle to matrix
%   rvec : axis-angle representation scaled by angle

theta = norm(rvec, 2);
k = rvec / theta;

k_cross = [0, -k(3), k(2);
    k(3), 0, -k(1);
    -k(2), k(1), 0];

I = eye(3);

rotation_matrix = I + sin(theta) * k_cross + (1 - cos(theta))*k_cross^2;

end

