function [ y ] = enhaceEdges( image, M, S )

ker_v = [1 2 1; 0 0 0; -1 -2 -1];
ker_h = ker_v';

sm = conv2(double(image), double(ker_v), 'same');
sn = conv2(double(image), double(ker_h), 'same');

grad_image = sqrt(sm.^2 + sn.^2);

filter = ones(M, M) ./ (M*M);
filtered_image = conv2(double(image), double(filter), 'same');
sa = conv2(double(grad_image), double(filter), 'same');

w = exp( -sa ./ S);

y = uint8((1 - w) .* double(image)  + w .* double(filtered_image));

end

