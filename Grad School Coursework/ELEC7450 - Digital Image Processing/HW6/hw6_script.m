% HW6
% Matt Boler

clc; clear all; close all
set(0, 'DefaultFigureColormap', gray(256))
%% Exercise 1

image = rgb2gray(imread('samford.png'));
ker = fspecial('disk', 6);

%image_F = fft2(im2double(image));
image_F = fft2(double(image));
ker_F = fft2(ker, size(image, 1), size(image, 2));
blurred_image = real(ifft2(ker_F .* image_F));
figure(1)
imagesc(blurred_image)
title("Image blurred via FT");



%% Exercise 2

% oneliner: ker, blurred_image -> inverse filtered image

FFT_inverted = fft2(blurred_image) ./ fft2(ker, size(blurred_image, 1), size(blurred_image, 2));
inverted_image = real(ifft2(FFT_inverted));
figure(2)
imagesc(inverted_image)
title("Inverse Filtered Image")

%% Exercise 3
sigma = 5;
noise = sigma * randn(size(image));
% There has to be a better way than this...
noise_image = blurred_image + noise;
figure(4)
imagesc(noise_image)
title("Image with Noise")

%% Exercise 4