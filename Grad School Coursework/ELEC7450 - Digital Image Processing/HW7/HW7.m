clc; clear all; close all;
set(0, 'DefaultFigureColormap', gray(256))

image = imread('samford.png');
imageGray = double(rgb2gray(image));

baseHistogram = hist(reshape(imageGray, 1, []), [-255:255]);
numPixels = sum(baseHistogram);
pBase = baseHistogram ./ numPixels;
lBase = log2(pBase + eps);
entropyBase = -sum(pBase .* lBase)

split1 = [0.5; 0.5];
[errorImage, predictedImage] = predictImage(imageGray, split1);

figure(1)
imagesc(imageGray);
title("Original Image");
figure(2)
imagesc(errorImage);
title("Prediction Error Image");

e = reshape(errorImage, 1, []);
errorHistogram = hist(e, [-255:255]);

p1 = errorHistogram ./ numPixels;
l1 = log2(p1 + eps);
entropy1 = -sum(p1 .* l1)

figure(3)
plot(-255:255, errorHistogram);
title("Error Image Histogram");

split2 = [0.45; 0.55];
[errorImage2, ~] = predictImage(imageGray, split2);
e2 = reshape(errorImage2, 1, []);
hist2 = hist(e2, [-255:255]);
p2 = hist2 ./ numPixels;
l2 = log2(p2 + eps);
entropy2 = -sum(p2 .* l2)

split3 = [0.55; 0.45];
[errorImage3, ~] = predictImage(imageGray, split3);
e3 = reshape(errorImage3, 1, []);
hist3 = hist(e3, [-255:255]);
p3 = hist3 ./ numPixels;
l3 = log2(p3 + eps);
entropy3 = -sum(p3 .* l3)

%% Jpeg Compression

for i = 0:5:100
   filename = ['images/image_', num2str(i), '.jpg'];
   imwrite(uint8(imageGray), filename, 'jpg', 'Quality', i);
end