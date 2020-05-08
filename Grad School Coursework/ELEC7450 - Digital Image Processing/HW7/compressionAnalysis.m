quality = 5:5:100;
size = [16.0, 24.9, 32.6, 39.2, 45.4, 50.9, 56.1, 60.4, 65.0, 69.0, 73.1, 78.1, 84.3, 91.6, 99.7, 111.9, 128.2, 155.4, 210.1, 342.3];

figure(1)
plot(quality, size)
xlabel("Quality Level")
ylabel("Image Size (kB)")
title("Image Quality vs File Size");

image = rgb2gray(imread('samford.png'));
mse = [];

for i = 5:5:100
   filename = ['images/image_', num2str(i), '.jpg'];
   compressedImage = imread(filename);
   mse(end+1) = mean( (image(:) - compressedImage(:)).^2 );
end

figure(2)
plot(quality, mse);
xlabel("Quality Level")
ylabel("Mean Squared Error")
title("Image Quality vs Mean Squared Error");

figure(3)
image20 = imread('images/image_20.jpg');
imshow(image20)
title("Image at Quality = 20");

figure(4)
image65 = imread('images/image_65.jpg');
imshow(image65)
title("Image at Quality = 65");