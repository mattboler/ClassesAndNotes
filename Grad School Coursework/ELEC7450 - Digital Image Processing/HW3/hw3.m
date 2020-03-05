clc; clear all; close all;

set(0, 'DefaultFigureColormap', gray(256))


%{
Matt Boler
HW3
%}

%% Exercise 1
%{
1. Form a convolution kernel with ker1 = ones(13,3)/39;. Take the 2-D FFT via
fft2(ker1,256,256);.
(a) Display the magnitude of the 2-D FFT. (Use abs.) You may want to use fftshift
to shift the origin in the Fourier domain to the center of the plot. Include the Fourier
magnitude image in your report.
(b) Give an explanation for why the magnitude plot appears as it does (in both coordinates)
in terms of simple Fourier transform pairs and properties.
(c) Use the FFT plot to explain what would happen to an image if it were convolved with
this 2-D signal ker1.
%}

%a
ker1 = ones(13,3)/39;
%ker1 = ones(25, 1)/25;
f_ker1 = fft2(ker1, 256, 256);
figure(1);
imagesc(fftshift(abs(f_ker1)));
title("FT of ker1");

%b

%c

%% Exercise 2
%{
2. Read in camera.tif and convert the integer values to a floating-point image with double.
(a) Take the 2-D FFT, and plot the magnitude. What does the plot look like? (Look very
closely. It isn’t all black.)
(b) Use the dynamic range compression formula (also used in Project 2) to transform the
magnitude coefficients. Scale so that the maximum value is 255 (or use imagesc or
imshow(I,[])), and plot. Describe the resulting plot.
(c) What strong features in the original image account for the most noticeable features in the
Fourier transform? (Think in terms of simple Fourier transform pairs and properties.)
%}

img = double(imread('camera.tif'));

%a
f_img = fft2(img);
figure(2)
imagesc(fftshift(abs(f_img)));
title("FT of Camera");

%b
figure(3)
imagesc(log(abs(fftshift(f_img))+10));
title("Log-Magnitude of FT of Camera");

%c

%% Exercise 3
%{
3. Multiply the FFT coefficients of the image and the kernel, and take the inverse FFT (ifft2).
(You may need to take the real part because the finite precision creates a very small imaginary
component.)
(a) Display the result, and describe it.
(b) What operation does this represent on the original cameraman image?
(c) Compare this to a linear convolution of the image with the kernel using conv2.
Include these images in your report.

%}

f_comb = f_ker1 .* f_img;
i_comb = real(ifft2(f_comb));

%a
figure(4)
imagesc(i_comb);
title("Inverse FT of ker1 * Camera");

%b

%c
comb = conv2(ker1, img);
figure(5);
imagesc(comb);
title("Convolved Filter and Image");

%% Exercise 4
img = zeros(1024, 1024);
tic
for i = 1:100
    img = rand(1024, 1024);
    fft2(img);
end
toc