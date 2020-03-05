ColorImage = imread('img_818.png'); %Read color image at t0
GrayImage = rgb2gray(ColorImage); %Convert to gray scale image

%SURF interest points
points_SURF = detectSURFFeatures(GrayImage);
subplot(221); imshow(GrayImage); hold on;
plot(points_SURF);
title('SURF interest points')

%BRISK interest points
points_BRISK = detectBRISKFeatures(GrayImage);
subplot(222); imshow(GrayImage); hold on;
plot(points_BRISK);
title('BRISK interest points')

%FAST interest points
points_FAST = detectFASTFeatures(GrayImage);
subplot(223); imshow(GrayImage); hold on;
plot(points_FAST);
title('FAST interest points')

%Harris interest points
points_Harris = detectHarrisFeatures(GrayImage);
subplot(224); imshow(GrayImage); hold on;
plot(points_Harris);
title('Harris interest points')