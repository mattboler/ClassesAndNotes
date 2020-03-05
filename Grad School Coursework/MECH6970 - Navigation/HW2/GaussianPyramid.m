ColorImage = imread('Image_t0.jpg'); %Read color image at t0
I = rgb2gray(ColorImage); %Covert to gray scale image
Is = im2single(I); %Convert to single
IsPyr = Is; %1st level of Pyramid
pyrLevels =4; %Number of Pyramid levels
subplot(1,pyrLevels,1) 
imshow(Is) %Show 1st level
title(['GPyrLevel: ', num2str(1)])
for i=2:1:pyrLevels
    Is = imgaussfilt(Is,2^(i-1)); %Gaussian filtered with variance 2^(i-1)
    IsPyr = [IsPyr Is];
    subplot(1,pyrLevels,i)
    imshow(Is); %Show ith level
    title(['GPyrLevel: ',num2str(i)])
end
