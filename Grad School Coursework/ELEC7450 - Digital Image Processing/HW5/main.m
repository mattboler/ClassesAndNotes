image1 = imread('books1.jpg');
image2 = imread('books2.jpg');

%% Part 2: Apply Affine Transform
disp("Part 2:")
A = [.5 .5; .5 -1.5];
b = [0; 0];

transformedImage = affineTransformBackward(A, b, image1);
figure(1)
imshow(image1);

figure(2)
imshow(transformedImage);

%% Part 3: Find Transform Between Images
disp("Part 3:")
%{

Points to be found:
Yellow book:
1: Top left corner of P
2: Top right corner of rightmost N
3: Inner-left corner of leftmost T
4. Inner-right corner of rightmost T
5. Bottom corner of blurb

Black book:
1. Top left of capital E
2. Inner right corner of rightmost lowercase e
3. Top left of F in 'FOR'
4. Top left of left E in 'ENGINEERS'
5. Inner right corner of rightmost T in 'SCIENTISTS'
%}

[h1, w1, d1] = size(image1);
du1 = floor(w1/2);
dv1 = floor(h1/2);
u1 = [485, 756, 622, 696, 437, 114, 352, 99, 152, 240];
v1 = [90, 97, 105, 106, 270, 159, 238, 274, 287, 378];
% correct coordinates
u1c = u1 - du1;
v1c = v1 - dv1;

[h2, w2, d2] = size(image2);
du2 = floor(w2/2);
dv2 = floor(h2/2);
u2 = [464 712 581 651 336 124 271 59 89 110];
v2 = [123 165 156 166 240 100 190 167 184 262];
% correct coordinates
u2c = u2 - du2;
v2c = v2 - dv2;

% u2 = a*u1 + b*u2 + c*1
% v2 = d*u1 + e*u2 + f*1
X = [u1', v1', ones(size(u1))'];
Yu = u2';
Yv = v2';

T = [(X \ Yu)'; (X \ Yv)'; 0 0 1];

%% Part 4: Apply transform
disp("Part 4:")
simulatedImage2 = affineTransformBackward(T(1:2, 1:2), T(1:2, 3), image1);
figure(3)
imshow(image2)
title("Image 2")

figure(4)
imshow(simulatedImage2)
title("Simulated Image 2")