function [out] = error_ellipse(covariance, k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

theta = linspace(0, 2*pi, 100);
a = k * [cos(theta); sin(theta)];

[V, D] = eig(covariance);
A = D^(-1/2) * V;
A_inv = pinv(A);

out = A_inv * a;

end

