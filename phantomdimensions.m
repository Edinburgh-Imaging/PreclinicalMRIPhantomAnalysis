function [analyseddimmasks, locationsArray, xdimension, ydimension] =...
    phantomdimensions(imagename, imagefolder, vol, plane, hpixsize,...
    vpixsize, centralslice, ni, nf, scalefactor, thresholdf, sigma,...
    saveimages, method)
% Companion function of the PhantomAnalysisGUI. It estimates the internal
% dimensions of the LEGO phantom by measuring the distance between 
% opposing edges at 4 different locations per direction in the imaging 
% plane and per slice and selects the mode value for each dimension. There
% are 3 methods available for identifying the boundary of the phantom:
% a. Canny edge detection.
% b. k-means clustering with 2 clusters: one cluster represents the phantom
%    and the other the background noise.
% c. k-means clustering with 3 clusters: two clusters represent the phantom
%    and the other the background noise.
% All three methods are accompanied by other morphological operations.
%
% >> [analyseddimmasks, locationsArray, xdimension, ydimension] =...
%    phantomdimensions(imagename, imagefolder, vol, plane, hpixsize,...
%    vpixsize, centralslice, ni, nf, scalefactor, thresholdf, sigma,...
%    saveimages, method)
%
% Variable Dictionary:
% --------------------
% imagename          input  The name of the image.
% imagefolder        input  The folder where a new folder will be created 
%                           for saving the analysed images.
% vol                input  A 3D matrix version of the MRI sequence. This
%                           is created by PhantomAnalysisGUI.
% plane              input  The MRI imaging plane the image was acquired.
% hpixsize           input  The size of the voxel in X direction.
% vpixsize           input  The size of the voxel in Y direction.
% centralslice       input  The central slice of the 3D matrix.
% ni                 input  The first slice to analyse.
% nf                 input  The last slice to analyse.
% scalefactor        input  The factor to scale the images for subpixel
%                           analysis.
% thresholdf         input  Factor to multiply Canny detector's thresholds.
% sigma              input  Standard deviation of Canny detector's Gaussian
%                           blurring.
% saveimages         input  Equals 1 if user wants to save analysed images
%                           and 0 if not.
% method             input  The method for identifying the phantom's
%                           boundary: 'Canny', 'kmeans2' or 'kmeans3'.
% analyseddimmasks   output A 3D matrix of the analysed masks.
% locationsArray     output A matrix containing the locations where
%                           vectors and edges overlap.
% xdimension         output The estimated distance (in mm) at the 
%                           horizontal direction.
% ydimension         output The estimated distance (in mm) at the vertical
%                           direction.
%
% Last Modified: 21 July 2016
% Copyright (c) 2016, Xenios Milidonis

% tic % stopwatch measuring elapsed time starts here

% Create new folder to save analysed images, if required.
if saveimages == 1
    mkdir(imagefolder,...
        sprintf('%s_phantomdimensions', imagename(1:end - 4)));
    savefolder = [imagefolder,...
        sprintf('%s_phantomdimensions', imagename(1:end - 4))];
else % must also define an alternative for code to work
    savefolder = 'none';
end

% % Determine rotation or not.
% switch plane
%     case {'Axial', 'Coronal', 'Coronal90', 'Sagittal90'}
%         rotate = 'no';
%     case {'Axial90', 'Sagittal'}
%         rotate = 'yes';
% end

% Get the central slice.
central = vol(:, :, centralslice);

% Call kmeanssegmentation.m to get the binarised central slice using 
% k-means clustering. Two white clusters out of a total of 3 is the best 
% option for phantom images with no uniform illumination.
bwcentral = kmeanssegmentation(central, 3, 2);

% If multiple objects are found, try to attach them.
bwcentral = attachobjects(bwcentral);

% Get the position of the centroid and round it up.
[centroidr, centroidc] = centroid(bwcentral);
centroidr = round(centroidr);
centroidc = round(centroidc);

% Estimate the horizontal dimension of the phantom using a row vector.
% Add 1 pixel to take into accound both pixels at the edges.
rowvector = bwcentral(centroidr, :);
[~, rowleftcol]  = find(rowvector, 1, 'first');
[~, rowrightcol] = find(rowvector, 1, 'last');
rowpixels = rowrightcol - rowleftcol + 1;

% Multiply by the scale factor.
horpixels = round(rowpixels * scalefactor);

% Call scalestack.m to to get the scaled images for sub-pixel analysis.
scaledseq = scalestack(vol, ni, nf, scalefactor, 'bicubic');

% Canny detector's sigma should be changed for sub-pixel analysis.
sigma = sigma * scalefactor / 2;

% Call savedimensions.m to save analysed images, if required.
if saveimages == 1
    savedimensions(savefolder, ni, nf, centroidr, centroidc, central,...
    scaledseq);
end

% Show in MATLAB's Command Window the headers for the measurements.
headerformatspec = '%9s %9s %9s %9s %9s %9s %9s %9s %9s \n';
fprintf(headerformatspec, 'Slice',...
    'HorU1(mm)', 'HorU2(mm)', 'HorD1(mm)', 'HorD2(mm)',...
    'VerL1(mm)', 'VerL2(mm)', 'VerR1(mm)', 'VerR2(mm)');

% Create 8 column arrays to store the measurements.
[horuArray, hordArray, verlArray, verrArray, horuArray2, hordArray2,... 
    verlArray2, verrArray2] = deal(zeros(nf - ni + 1, 1));

% Create a matrix to store the positions of the 16 points per slice that
% will be used for visualisation in the GUI's 'Output Data' figure.
locationsArray = NaN(16, 2, nf - ni + 1);

% Create a marker image with a white border of 2 pixels width around the 
% image to be used for reconstruction.
% For k-means check if images have a strong bias field by comparing the 
% mean of the intensity in the upper and lower half of the central slice.
% if strcmp(rotate, 'no')
    markerborder = ones(size(central));
    markerborder(3:size(central, 1) - 2, 3:size(central, 2) - 2) = 0;  

    meanAbove =...
        mean(mean(central(centroidr - (size(central, 1)/5):centroidr, :)));
    meanBelow =...
        mean(mean(central(centroidr:centroidr + (size(central, 1)/5), :)));     
% else   % image size and centroid might be different
%     rotated = imrotate(central, -90);
%     [centroidr, centroidc] = centroid(rotated);
%     centroidr = round(centroidr);
%     centroidc = round(centroidc);
%     
%     markerborder = ones(size(rotated));
%     markerborder(3:size(rotated, 1) - 2, 3:size(rotated, 2) - 2) = 0; 
% 
%     meanAbove =...
%         mean(mean(rotated(centroidr - 50:centroidr, :)));
%     meanBelow =...
%         mean(mean(rotated(centroidr:centroidr + 50, :)));   
% end

% Resize markerborder for sub-pixel analysis with nearest-neighbor 
% interpolation to preserve pixels with value of 1.
markerborder = imresize(markerborder, scalefactor, 'nearest');

% Main loop for measuring the dimensions of the phantom.
for i = 1:(nf - ni + 1) % slice 1 is slice ni of original sequence                                            
    rotated = scaledseq(:, :, i);

%     % Rotate images by 90 degrees if required. For this case the centroid
%     % was re-calculated earlier.
%     if strcmp(rotate, 'yes')
%         rotated = imrotate(scaled, -90);
%     else
%         rotated = scaled;
%     end

    % Pre-processing steps for k-means.
    if strcmp(method, 'kmeans2') || strcmp(method, 'kmeans3')
        % Take the sqrt of the image if there is a strong bias field, i.e. 
        % the relative difference in the intensity of the two halves of the
        % phantom is larger than 35%.
        if (abs(meanAbove - meanBelow) / mean([meanAbove, meanBelow]))...
                > 0.35      
            rotated = sqrt(rotated);
        end
        
        % Blurring helps reduce noise. This is a default step in Canny.
        h = fspecial('gaussian', [3 3], 1);
        rotated = imfilter(rotated, h);       
    end
     
    % Apply the appropriate method to identify the boundary of the phantom.
    switch method
        case 'Canny'
            % Apply edge detection.
            [~, t] = edge(rotated, 'canny');
            [edged, ~] = edge(rotated, 'canny', t * thresholdf, sigma); 
            
        case 'kmeans2'
            % Call kmeanssegmentation.m to apply k-means clustering.
            segmented = kmeanssegmentation(rotated, 2, 1);
    
        case 'kmeans3'
            % Call kmeanssegmentation.m to apply k-means clustering.
            segmented = kmeanssegmentation(rotated, 3, 2);
            
            % Rarely, clusters 1 and 2 are exclusively background noise,
            % thus segmentation fails (selection includes background). In 
            % this case use k-means with k = 2.
            shapef = shapefactor(segmented, 'c');
            if shapef < 0.1
                segmented = kmeanssegmentation(rotated, 2, 1);
            end
    end
        
    % Post-processing steps for all methods.
    switch method
        case 'Canny'
            % Remove some of the unwanted pixels at the border of the image
            % and around the phantom.
            edged = +edged;
            boundarypixels = imreconstruct(markerborder, edged);
            cleaned = logical(edged - boundarypixels);
            cleaned = bwareaopen(cleaned, ceil(size(cleaned, 1) / 7));
            
        case {'kmeans2', 'kmeans3'}
            % Get the boundary of the segmented phantom. This is needed for
            % saving images of the analysis process only.
            edged = bwperim(segmented, 4);

            % Remove some of the unwanted pixels at the border of the image
            % and around the phantom.
            segmented = +segmented;
            boundarypixels = imreconstruct(markerborder, segmented);
            cleaned = logical(segmented - boundarypixels);
            cleaned = bwareaopen(cleaned, ceil(size(cleaned, 1) / 7));

            % Smooth the sides of the volumes using structural elements.
            % Canny does not require further smoothing.
            newcentroidr = round(centroidr * scalefactor);
            newcentroidc = round(centroidc * scalefactor);
            cleaned = fillandsmooth(cleaned, newcentroidr, newcentroidc,...
                'centroid');

            % Get the boundary of the cleaned phantom. This is the one
            % needed for calculating the dimensions.
            cleaned = bwperim(cleaned, 4);
    end
     
    % If sub-pixel analysis was used, resize the analysed slices back to
    % the original size to improve the speed of visualisation. To preserve
    % the detected boundaries some analysis is required.
    if scalefactor ~= 1
        % Use a retangular structural element to dilate the image.
        se = strel('rectangle', [ceil(scalefactor), ceil(scalefactor)]);
        analyseddimslice = imdilate(cleaned, se);
        
        % Resize back to original size using with nearest-neighbor 
        % interpolation to preserve pixels with value of 1.
        analyseddimslice = imresize(analyseddimslice, 1 / scalefactor,...
            'nearest');
        
        % Perform morphological thinning.
        analyseddimslice = bwmorph(analyseddimslice, 'thin', Inf);
        
        % Convert to double.
        analyseddimslice = +analyseddimslice;
    else
        % Convert to double.
        analyseddimslice = +cleaned;
    end
        
    % Create a 3D matrix version of the analysed slices. This will be
    % passed as a function output to be shown in the GUI's 'Output Data'
    % figure. 
    if i == 1                                       
        analyseddimmasks = analyseddimslice;
    else
        analyseddimmasks = cat(3, analyseddimmasks, analyseddimslice);
    end

    % Call vectorsdimensions.m to measure the dimensions at multiple 
    % locations. 
    [horuArray, hordArray, verlArray, verrArray, horuArray2,...
    hordArray2, verlArray2, verrArray2, locationsArray, shiftUD,...
    shiftUD2, shiftLR, shiftLR2] =...
    vectorsdimensions(cleaned, centroidr, centroidc, scalefactor,...
    horuArray, hordArray, verlArray, verrArray, horuArray2, hordArray2,...
    verlArray2, verrArray2, locationsArray, plane, hpixsize, vpixsize,...
    horpixels, i, method);

    % Call savedimensionsloop.m to save analysed slices, if required.
    if saveimages == 1
        savedimensionsloop(savefolder, ni, nf, centroidr, centroidc,...
        scalefactor, rotated, edged, cleaned, shiftUD, shiftUD2,...
        shiftLR, shiftLR2, i);
    end              
    
    % Show in MATLAB's Command Window the measurements.
    formatspec = '%9d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n';
    fprintf(formatspec, i + ni - 1,...
        horuArray(i), horuArray2(i), hordArray(i), hordArray2(i),...
        verlArray(i), verlArray2(i), verrArray(i), verrArray2(i));
end

% Calculate the mean, mode and median dimensions.
horArray = cat(1, horuArray, hordArray, horuArray2, hordArray2);
horArray = horArray(horArray ~= 0);    % remove zeros
horMedian = median(horArray);
horlargeoutliers = horArray > 1.1*horMedian;
horArray(horlargeoutliers) = [];       % remove large outliers
horsmalloutliers = horArray < 0.9*horMedian;
horArray(horsmalloutliers) = [];       % remove small outliers
horMean = mean(horArray);
horMedi = median(horArray);
horModeSmall = mode(horArray);         % find the largest mode:
horArrayWithoutSmallMode = horArray(horArray ~= horModeSmall);   
horModeLarge = mode(horArrayWithoutSmallMode);
if sum(horArray(:) == horModeSmall)...
        == sum(horArrayWithoutSmallMode(:) == horModeLarge)
    horMode = horModeLarge;
else
    horMode = horModeSmall;
end

verArray = cat(1, verlArray, verrArray, verlArray2, verrArray2);
verArray = verArray(verArray ~= 0);    % remove zeros
verMedian = median(verArray);
verlargeoutliers = verArray > 1.1*verMedian;
verArray(verlargeoutliers) = [];       % remove large outliers
versmalloutliers = verArray < 0.9*verMedian;
verArray(versmalloutliers) = [];       % remove small outliers
verMean = mean(verArray);
verMedi = median(verArray);
verModeSmall = mode(verArray);         % find the largest mode:
verArrayWithoutSmallMode = verArray(verArray ~= verModeSmall);   
verModeLarge = mode(verArrayWithoutSmallMode);
if sum(verArray(:) == verModeSmall)...
        == sum(verArrayWithoutSmallMode(:) == verModeLarge)
    verMode = verModeLarge;
else
    verMode = verModeSmall;
end

% Show in MATLAB's Command Window the final dimensions.
fprintf(['Horizontal: mean = %7.3f mm, mode = %7.3f mm, '...
    'median = %7.3f mm\n'], horMean, horMode, horMedi);
fprintf(['  Vertical: mean = %7.3f mm, mode = %7.3f mm, '...
    'median = %7.3f mm\n'], verMean, verMode, verMedi);

% Pass this as an output to be shown in the GUI's text boxes.
xdimension = horMode;
ydimension = verMode;

% toc % stopwatch measuring elapsed time stops here

