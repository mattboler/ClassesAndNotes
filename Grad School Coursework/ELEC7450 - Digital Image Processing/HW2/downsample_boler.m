function [downsampled] = downsample_boler(image)
% My plan: 4 separate halved images, add and average
[height, width] = size(image);
u_idx = 1:2:width;
v_idx = 1:2:height;
img_topleft = image(u_idx, v_idx);
img_topright = image(u_idx+1, v_idx);
img_botleft = image(u_idx, v_idx+1);
img_botright = image(u_idx+1, v_idx+1);
downsampled = (img_topleft + img_topright + img_botleft + img_botright) ./ 4;
end

