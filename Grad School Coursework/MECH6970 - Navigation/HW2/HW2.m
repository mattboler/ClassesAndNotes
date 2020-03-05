clc; clear all; close all;

%{
Matt Boler
Fundamentals of Nav and Guidance
HW2
%}

%% Problem 1
p1_image = rgb2gray(imread('test-image.jpg'));
thresholds = linspace(500, 5000);
detections = zeros(size(thresholds));

for i = 1:size(thresholds,2)
    I = detectSURFFeatures(p1_image, 'MetricThreshold', thresholds(i));
    detections(1,i) = size(I,1);
    clear I;
end

figure(1)
plot(thresholds, detections);
title("Number of Features Detected vs Metric Threshold");
xlabel("Metric Threshold");
ylabel("Features Detected");

%% Problem 2
p2_im1 = rgb2gray(imread('Image_t0.jpg'));
p2_im2 = rgb2gray(imread('Image_t1.jpg'));

pts1 = detectSURFFeatures(p2_im1);
[f1, vpts1] = extractFeatures(p2_im1, pts1);

pts2 = detectSURFFeatures(p2_im2);
[f2, vpts2] = extractFeatures(p2_im2, pts2);

thresholds = linspace(1, 50, 1000);
ratios = linspace(0.01,1);

matches_1 = zeros(size(thresholds));

for i = 1:size(thresholds, 2)
    idxPairs = matchFeatures(f1, f2, 'Unique', true, 'MatchThreshold', thresholds(i));
    matches_1(1,i) = size(vpts1(idxPairs(:,1)), 1);
end

figure(2)
plot(thresholds, matches_1);
title("Number of Features Matched vs Match Threshold");
xlabel("Match Threshold");
ylabel("Features Matched");

matches_2 = zeros(size(matches_1));

for i = 1:size(ratios, 2)
   idxPairs = matchFeatures(f1, f2, 'Unique', true, 'MaxRatio', ratios(i));
   matches_2(1,i) = size(vpts1(idxPairs(:,1)),1);
end

figure(3)
plot(thresholds, matches_2);
title("Number of Features Matched vs Max Ratio");
xlabel("Max Ratio");
ylabel("Features Matched");

%% Problem 3

p3_im1 = rgb2gray(imread('reference_map.jpg'));
p3_im2 = rgb2gray(imread('img_132.png'));

max_features_matched = -1;
best_det_thresh = -1;
best_match_thresh = -1;
best_match_ratio = -1;

det_thresholds = linspace(500, 2000, 10);
match_thresholds = linspace(1, 10, 10);
match_ratios = linspace(0.01,1, 10);

for i = 1:size(det_thresholds, 2)
    det_thresh = det_thresholds(i);
    pts1 = detectSURFFeatures(p3_im1, 'MetricThreshold', det_thresh);
    [f1, vpts1] = extractFeatures(p3_im1, pts1);
    
    pts2 = detectSURFFeatures(p3_im2, 'MetricThreshold', det_thresh);
    [f2, vpts2] = extractFeatures(p3_im2, pts2);
    
    for j = 1:size(match_thresholds, 2)
        match_thresh = match_thresholds(j);
        for k = 1:size(match_ratios, 2)
            match_ratio = match_ratios(k);
            
            % Test parameters
            idxPairs = matchFeatures(f1, f2, 'Unique', true, 'MatchThreshold', match_thresh, 'MaxRatio', match_ratio);
            temp_matches = size(vpts1(idxPairs(:,1)), 1);
            
            if temp_matches > max_features_matched
               max_features_matched  = temp_matches;
               best_det_thresh = det_thresh;
               best_match_thresh = match_thresh;
               best_match_ratio = match_ratio;
            end
            
        end
    end
end

disp("Best Detection Threshold: ")
disp(best_det_thresh)
disp("Best Matching Threshold: ")
disp(best_match_thresh)
disp("Best Matching Ratio: ")
disp(best_match_ratio)
disp("Features Matched: ")
disp(max_features_matched)

%% Bonus

p3_im1 = rgb2gray(imread('Image_t0.jpg'));
p3_im2 = rgb2gray(imread('Image_t1.jpg'));

max_features_matched = -1;
best_det_thresh = -1;
best_match_thresh = -1;
best_match_ratio = -1;

det_thresholds = linspace(500, 2000, 10);
match_thresholds = linspace(1, 10, 10);
match_ratios = linspace(0.01,1, 10);

for i = 1:size(det_thresholds, 2)
    det_thresh = det_thresholds(i);
    pts1 = detectSURFFeatures(p3_im1, 'MetricThreshold', det_thresh);
    [f1, vpts1] = extractFeatures(p3_im1, pts1);
    
    pts2 = detectSURFFeatures(p3_im2, 'MetricThreshold', det_thresh);
    [f2, vpts2] = extractFeatures(p3_im2, pts2);
    
    for j = 1:size(match_thresholds, 2)
        match_thresh = match_thresholds(j);
        for k = 1:size(match_ratios, 2)
            match_ratio = match_ratios(k);
            
            % Test parameters
            [idxPairs, metric] = matchFeatures(f1, f2, 'Unique', true, 'MatchThreshold', match_thresh, 'MaxRatio', match_ratio);
            temp_matches = size(vpts1(idxPairs(:,1)), 1);
            temp_avg_metric = sum(metric) / temp_matches;
            
            if temp_matches > max_features_matched
               max_features_matched  = temp_matches;
               best_det_thresh = det_thresh;
               best_match_thresh = match_thresh;
               best_match_ratio = match_ratio;
            end
            
        end
    end
end

disp("BONUS")

disp("Best Detection Threshold: ")
disp(best_det_thresh)
disp("Best Matching Threshold: ")
disp(best_match_thresh)
disp("Best Matching Ratio: ")
disp(best_match_ratio)
disp("Features Matched: ")
disp(max_features_matched)


