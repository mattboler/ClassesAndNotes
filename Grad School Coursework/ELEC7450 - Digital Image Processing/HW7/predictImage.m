function [errorImage, predictedImage] = predictImage(image, split)
%   split is [weight1, weight2]

predictedImage = zeros(size(image), 'double');
a = split(1);
b = split(2);

% NOTE: Function is undefined for m = 1 or n = 1 b/c (out of range)

for m = 2:size(image, 1)
    for n = 2:size(image,2)
        predictedImage(m, n) = a * image(m, n-1) + b * image(m-1, n);
    end
end

errorImage = image - predictedImage;

end

