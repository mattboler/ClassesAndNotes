clear all
close all

%% ---- PARAMETERS ----
enableWarp = 0; % 0: Disables Warping, 1: Enables Warping
numberFeaturesDetected = 1000; % Number of points to select with SURF algorithm
enableSelectingNewLocation = 0; % 0: Uses default crop location, 1: Allows user to select new location with rectange
enableCrop = 1; % 0: use raw image, 1: crop raw image
enableCalculatedCamera = 1; % 0: disable calculated camera, 1: enable calculated camera

%% ---- STORAGE ----
num_features = [200, 400, 600, 800, 1000];
calculated_orientations = {};
simulated_orientations = {};
calculated_locations = {};
simulated_locations = {};

for index = 1:length(num_features)
    numberFeaturesDetected = num_features(index);
    
    %% ---- INTEREST POINT COMPUTATION FOR EACH IMAGE ----
    ColorImage_t0 = imread('Google-maps-ref.jpg'); %Read color image at t0
    I_t0 = rgb2gray(ColorImage_t0); %Convert to gray scale image

    points_t0 = selectStrongest(detectSURFFeatures(I_t0),numberFeaturesDetected); %Detect SURF features at t0

    if enableSelectingNewLocation == 1 
        figure, imshow(ColorImage_t0)
        map_sub_selection = getrect()
        print("Click and drag in figure to select image crop")
    else
        map_sub_selection = 1.0e+03 *[1.0198    0.6388    0.3345    0.3165];
    end

    ColorImage_t1 = imread('Google-maps-ref.jpg'); %Rinse and repeat for image at t1
    if enableCrop == 1
        ColorImage_t1 = imcrop(ColorImage_t1,map_sub_selection);
    end
    I_t1 = rgb2gray(ColorImage_t1);

    % Perform perspective warp on cropped image
    if enableWarp == 1
        theta = 10;
        scale = 0.0001
        tm = [cosd(theta) -sind(theta) scale; ...
            sind(theta) cosd(theta) 0; ...
            scale 0 1];
        tform = projective2d(tm);
        I_t1 = imwarp(I_t1,tform);
        figure(5);imshow(I_t1);
    end

    % points_t1 = detectSURFFeatures(I_t1);

    points_t1 = selectStrongest(detectSURFFeatures(I_t1),numberFeaturesDetected);

    %EXTRACT FEATURES FROM INTEREST POINTS 
    [features_t0,validPoints_t0] = extractFeatures(I_t0,points_t0); %Extract features from interest point at t0
    [features_t1,validPoints_t1] = extractFeatures(I_t1,points_t1); %Same for t1

    %MATCH EXTRACTED FEATURES 
    indexPairs = matchFeatures(features_t0, features_t1); %Match interest point features at t0 with t1

    %DISPLAY MATCHED INTEREST POINTS
    matchedPoints_t0 = validPoints_t0(indexPairs(:, 1), :); %Retrieve the locations of matched points at t0
    %subplot(131); imshow(I_t0); hold on; %Display image at t0
    %plot(matchedPoints_t0); %Overlay matched points at t0

    matchedPoints_t1 = validPoints_t1(indexPairs(:, 2), :); %Same for t1 
    %subplot(132); imshow(I_t1); hold on;
    %plot(matchedPoints_t1);

    % [fLMedS, inliers] = estimateFundamentalMatrix(matchedPoints_t0,matchedPoints_t1,'NumTrials',5000);

    % 
    % figure(2); ax = axes;
    % showMatchedFeatures(I_t0, I_t1, matchedPoints_t0(inliers,:), matchedPoints_t1(inliers,:), 'falsecolor','Parent',ax); %Side-by-side of both images with matched points
    % title(ax, 'Candidate point matches');
    % legend(ax, 'Matched points @t0','Matched points @t1');

    P_t0 = double(matchedPoints_t0.Location)'; %Matched points from image at t0
    zP_t0 = ones(1,size(P_t0,2));
    P_t0= [P_t0; zP_t0]'; %Make each matched point a 1x3 vector 

    P_t1 = double(matchedPoints_t1.Location)'; %Repeat for image at t1
    zP_t1 = ones(1,size(P_t1,2));
    P_t1= [P_t1; zP_t1]';

    z = zeros(length(P_t0(:,1)),1);
    F = [P_t0(:,1),P_t0(:,2),P_t0(:,3),z,z,z,-P_t0(:,1).*P_t1(:,1),-P_t0(:,2).*P_t1(:,1);
        z,z,z,P_t0(:,1),P_t0(:,2),P_t0(:,3),-P_t0(:,1).*P_t1(:,2),-P_t0(:,2).*P_t1(:,2)]; %Compute F Matrix 

    B = [P_t1(:,1);
        P_t1(:,2)]; %Compute B matrix

    Gtilde = F \ B; %Solve for Gtilde matrix
    Ghat(1,:)=Gtilde(1:3); 
    Ghat(2,:)=Gtilde(4:6);
    Ghat(3,1:2)=Gtilde(7:8); Ghat(3,3)=1; %Reshape Gtilde into Ghat

    [U,L,V] = svd(Ghat); %SVD of Ghat

    if(L(1,1) > L(2,2) && L(2,2) > L(3,3))
        d = L(2,2); %Distance of plane
        %Computation of nhat using constraint that n'.P_t0(i) > 0, using i=3
        ntilde=[sqrt((L(1,1)^2-L(2,2)^2)/(L(1,1)^2-L(3,3)^2));
            0;
            sqrt((L(2,2)^2-L(3,3)^2)/(L(1,1)^2-L(3,3)^2))];
        shat= zeros(1,4);
        if(ntilde'*[1 0 0;0 1 0;0 0 1]*V'*P_t0(3,:)' > 0)
            nhat = [1 0 0;0 1 0;0 0 1]*ntilde;
            end
        if(ntilde'*[-1 0 0;0 1 0;0 0 -1]*V'*P_t0(3,:)' > 0)
            nhat = [-1 0 0;0 1 0;0 0 -1]*ntilde;
            end
        if(ntilde'*[-1 0 0;0 1 0;0 0 1]*V'*P_t0(3,:)' > 0)
            nhat = [-1 0 0;0 1 0;0 0 1]*ntilde;
            end
        if(ntilde'*[1 0 0;0 1 0;0 0 -1]*V'*P_t0(3,:)' > 0)
            nhat = [1 0 0;0 1 0;0 0 -1]*ntilde;
        end
        %Computation of xhat based on nhat
        xhat=(L(1,1)-L(3,3))*[nhat(1);0;-nhat(3)];
    end
    %Computation of Rhat
    Rhat=zeros(3,3);
    Rhat(2,:)=[0 1 0];
    Rhat(:,2)=Rhat(2,:)';
    Rhat(1,1)=(L(2,2)^2+L(1,1)*L(3,3))/((L(1,1)+L(3,3))*L(2,2));
    Rhat(3,3)=Rhat(1,1);
    Rhat(3,1)=sign(nhat(1))*sign(nhat(3))*sqrt((L(1,1)^2-L(2,2)^2)*(L(2,2)^2-L(3,3)^2))/((L(1,1)+L(3,3))*L(2,2));
    Rhat(1,3)=-Rhat(3,1);

    n = V'*nhat; %See slides
    x01 = U*xhat;
    R01 = U*Rhat*V';

    %% ---- Estimate Position Based on Computer Vision Toolbox ----
    focalLength = 400; % arbitrary value, only adjusts z value in orthogonal case

    % Determine Camera Paramters
    [rows1, columns1] = size(I_t1);
    IntrinsicMatrix = [focalLength,0,0;0,focalLength,0;columns1/2,rows1/2,1];
    radialDistortion = [0,0,0]; 
    cameraParams = cameraParameters('IntrinsicMatrix',IntrinsicMatrix,'RadialDistortion',radialDistortion); 

    % Reformat calculated values
    calculatedLocation = transpose(x01);
    calculatedLocation(3) = calculatedLocation(3)*focalLength*-1;
    calculatedOrientation = R01;

    % Determine World Orientation and Location
    worldPoints = [matchedPoints_t0.Location zeros(length(matchedPoints_t0.Location),1)];
    imagePoints = matchedPoints_t1.Location;
    [worldOrientation,worldLocation] = estimateWorldCameraPose(imagePoints ,worldPoints,cameraParams);


    % plot arbitrary world camera
    [rows0, columns0] = size(I_t0);
    worldCameraHeight = -focalLength; % arbitrary values

    % plot calculated camera
    calculatedWorldLocation = [calculatedLocation(1)+columns0/2 calculatedLocation(2)+rows0/2 calculatedLocation(3)];

    %% ---- Compare results of two methods ----
    simulatedLocation = [worldLocation(1)-columns0/2 worldLocation(2)-rows0/2 worldLocation(3)];
    calculatedLocation;

    simulatedOrientation = worldOrientation;
    calculatedOrientation;
    
    calculated_orientations{index} = calculatedOrientation;
    simulated_orientations{index} = simulatedOrientation;
    calculated_locations{index} = calculatedLocation;
    simulated_locations{index} = simulatedLocation;

end