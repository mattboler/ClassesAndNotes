%{
Project 2

This project focuses on sampling 2D signals.
%}

clc; clear all; close all;

%% Exercise 1
image_1 = zeros(1024, 1024);
image_1(500:502, :) = 255;

image_z = image_1(256:(256+511), 256:(256+511));

figure(1);
imagesc(image_z);
colormap(gray(256));
truesize
title("Horizontal Line")

image_rotated = imrotate(image_1, 25, 'bicubic');
image_rz = image_rotated(400:(400+511), 256:(256+511));
figure(2);
imagesc(image_rz);
colormap(gray(256));
truesize
title("Horizontal Line Rotated")

%{
Commentary: 
%}

%% Exercise 2

figure(3)
%imagesc(repmat(log(abs(fftshift(fft2(image_z)))+10.01), 1, 1))
imagesc(repmat(abs(fftshift(fft2(image_z))), 3, 3))
truesize
colormap(gray(512))
title("Fourier Magnitude of Unrotated Image");

figure(4)
%imagesc(repmat(abs(fftshift(fft2(image_zr)))+10.01), 1, 1))
imagesc(repmat(abs(fftshift(fft2(image_rz))), 3, 3))
truesize
colormap(gray(512))
title("Fourier Magnitude of Rotated Image");

%% Exercise 3

idxs = 1:2:512;

image_z_aliased = image_z(idxs, idxs);
image_rz_aliased = image_rz(idxs, idxs);

figure(5)
imagesc(image_z_aliased)
colormap(gray(256))
truesize
title("Aliased Nonrotated Image");

figure(6)
imagesc(image_rz_aliased)
colormap(gray(256))
truesize
title("Aliased Rotated image")

figure(7)
imagesc(repmat(abs(fftshift(fft2(image_z_aliased))), 3, 3))
truesize
colormap(gray(512))
title("Fourier Magnitude of Aliased, Unrotated Image");

figure(8)
imagesc(repmat(log(abs(fftshift(fft2(image_rz_aliased)))+10), 3, 3))
truesize
colormap(gray(512))
title("Fourier Magnitude of Aliased, Rotated Image");

%% Exercise 4

image_z_downsampled = downsample_boler(image_z);
image_rz_downsampled = downsample_boler(image_rz);

figure(9)
imagesc(image_rz_downsampled)
truesize
colormap(gray(512))
title("Downsampled Rotated Image");

figure(10)
imagesc(repmat(log(abs(fftshift(fft2(image_rz_downsampled)))+10), 3, 3))
truesize
colormap(gray(512))
title("Fourier Magnitude of Downsampled, Rotated Image");

%% Exercise 5
indm = 1:512;
indn = 1:512;

[x, y] = meshgrid(indm, indn);
m_factor = 0.65 * pi * y;
n_factor = 0.2 * pi * x;

f_image = cos(m_factor + n_factor);
% Obvious frequencies are wx = 0.65*pi rad/sample, wy = 0.2*pi rad/sample

idxs = 1:2:512;
f_image_aliased = f_image(idxs, idxs);

figure(11)
imagesc(f_image)
colormap(gray(512))
truesize
title("Sinusoid Image");

figure(12)
imagesc(f_image_aliased)
colormap(gray(512))
truesize
title("Aliased Sinusoid Image");

figure(13)
imagesc(repmat(abs(fftshift(fft2(f_image))), 3, 3))
colormap(gray(512))
truesize
title("Fourier Magnitude of Sinusoid Image");

%{
Commentary: 
peak at (423, 308) -> dm = 167, dn = 52
delta = 2*pi/512
freq_m = delta * dm = 0.65*pi
freq_n = delta*dn = 0.2*pi
%}

figure(14)
imagesc(repmat(abs(fftshift(fft2(f_image_aliased))), 3, 3))
colormap(gray(512))
truesize
title("Fourier Magnitude of Aliased Sinusoid Image");

%{
Commentary: 
peak at (39, 180) -> dm = -89, dn = 53

delta = 2*pi / 256
freq_m = delta * dm = -2.184 = -0.69*pi 
freq_n = delta * dn = 1.276 = 0.406*pi

hidden peak at (551, 436) -> dm = 167, dn = 52
freq_m = delta * dm = 1.30*pi = 2*original b/c scaling 
freq_n = dleta * dn = 0.406*pi = 2*original b/c scaling

-> IMPORTANT
peak_1_m is actually at 2*pi - hidden_freq_m!
peak_2_m is where expected because not large enough for aliasing to overlap
%}
