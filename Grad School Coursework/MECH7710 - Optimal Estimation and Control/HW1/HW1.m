clc; clear all; close all;

%{
MECH 7710 HW1
%}

%% 1: Generate discrete PDFs for the following:

f_die = (1/6)*[1, 1, 1, 1, 1, 1];
values_die = 1:6;
mean_die = mean(values_die);
var_die = sum((values_die - mean_die).^2) / length(values_die);

f_shifted_die = [0 0 0 f_die];
values_shifted_die = 4:9;
mean_shifted_die = mean(values_shifted_die);
var_shifted_die = sum((values_shifted_die - mean_shifted_die).^2) / length(values_shifted_die);

f_funny_die = (1/6)*[2, 0, 3, 0, 1, 0];
values_funny_die = [1 1 3 3 3 5];
mean_funny_die = mean(values_funny_die);
var_funny_die = sum((values_funny_die - mean_funny_die).^2) / length(values_funny_die);

figure(1)

% 6 dice numbered 1, 2, 3, 4, 5, 6
pdf_1a = conv(f_die, ...
    conv([0,f_die], ...
    conv([0, f_die], ...
    conv([0, f_die], ...
    conv([0, f_die], [0, f_die])))));
range_1a = 1:length(pdf_1a);
gauss_1a = normpdf(range_1a, 6*mean_die, sqrt(6*var_die));
subplot(2,2,1)
plot(pdf_1a);
hold on;
plot(gauss_1a);
title('1-A');
legend("PDF", "Normal");


% 6 dice numbered 4, 5, 6, 7, 8, 9
pdf_1b = conv(f_shifted_die, ...
    conv([0, f_shifted_die], ...
    conv([0, f_shifted_die], ...
    conv([0, f_shifted_die], ...
    conv([0, f_shifted_die], [0, f_shifted_die])))));
range_1b = 1:length(pdf_1b);
gauss_1b = normpdf(range_1b, 6*mean_shifted_die, sqrt(6*var_shifted_die));
subplot(2,2,2)
plot(pdf_1b)
hold on
plot(gauss_1b)
title('1-B')
legend("PDF", "Normal")

% 6 dice numbered 1 1 3 3 3 5
pdf_1c = conv(f_funny_die, ...
    conv([0,f_funny_die], ...
    conv([0, f_funny_die], ...
    conv([0, f_funny_die], ...
    conv([0, f_funny_die], [0, f_funny_die])))));
range_1c = 1:length(pdf_1c);
gauss_1c = normpdf(range_1c, 6*mean_funny_die, sqrt(6*var_funny_die));
subplot(2,2,3)
plot(pdf_1c)
hold on
plot(gauss_1c)
title("1-C");
legend("PDF", "Normal")

% 3 dice numbered 1 2 3 4 5 6 and 3 dice numbered 1 1 3 3 3 5
pdf_1d = conv(f_die, ...
    conv([0, f_die], ...
    conv([0, f_die], ...
    conv([0, f_funny_die], ...
    conv([0, f_funny_die], [0, f_funny_die])))));
range_1d = 1:length(pdf_1d);
gauss_1d = normpdf(range_1d, 3*mean_die + 3*mean_funny_die, sqrt(3*var_die + 3*mean_funny_die));
subplot(2,2,4)
plot(pdf_1d)
hold on
plot(gauss_1d)
title("1-D")
legend("PDF", "Normal")

%% 2
% 