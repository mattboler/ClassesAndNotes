%INTEREST POINT COMPUTATION FOR EACH IMAGE
ColorImage_t0 = imread('Image_t0.jpg'); %Read color image at t0
I_t0 = rgb2gray(ColorImage_t0); %Covert to gray scale image

points_t0 = detectSURFFeatures(I_t0); %Detect SURF features at t0
subplot(131); imshow(I_t0); hold on; %Display image at t0
plot(points_t0); %Overlay interest points at t0
 
ColorImage_t1 = imread('Image_t1.jpg'); %Rinse and repeat for image at t1
I_t1 = rgb2gray(ColorImage_t1);
points_t1 = detectSURFFeatures(I_t1);
subplot(132); imshow(I_t1); hold on;
plot(points_t1);

%EXTRACT FEATURES FROM INTEREST POINTS 
[features_t0,validPoints_t0] = extractFeatures(I_t0,points_t0); %Extract features from interest point at t0
[features_t1,validPoints_t1] = extractFeatures(I_t1,points_t1); %Same for t1

%MATCH EXTRACTED FEATURES 
indexPairs = matchFeatures(features_t0, features_t1); %Match interest point features at t0 with t1
 
%DISPLAY MATCHED INTEREST POINTS
matchedPoints_t0 = validPoints_t0(indexPairs(:, 1), :); %Retrieve the locations of matched points at t0
matchedPoints_t1 = validPoints_t1(indexPairs(:, 2), :); %Same for t1 
subplot(133); showMatchedFeatures(I_t0, I_t1, matchedPoints_t0, matchedPoints_t1); %Overlay of both images with matched points


 