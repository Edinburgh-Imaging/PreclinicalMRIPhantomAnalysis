function [analysedvolmasks, volume] = phantomvolume...
    (imagename, imagefolder, vol, plane, hpixsize, vpixsize,...
    slicethickness, centralslice, ni, nf, thresholdf, sigma,...
    saveimages, method)
% Companion function of the PhantomAnalysisGUI. It estimates the volume of
% the central frustum-shaped compartment of the LEGO phantom by identifying
% and then summing the corresponding region in each slice (Cavalieri's 
% principle). There are 3 methods available for segmenting the region of 
% interest on each slice:
% a. Canny edge detection.
% b. k-means clustering with 2 clusters: one cluster represents the phantom
%    and the other the background noise.
% c. k-means clustering with 3 clusters: two clusters represent the phantom
%    and the other the background noise.
% All three methods are accompanied by other morphological operations.
%
% >> [analysedvolmasks, volume] = phantomvolume...
%    (imagename, imagefolder, vol, plane, hpixsize, vpixsize,...
%    slicethickness, centralslice, ni, nf, thresholdf, sigma,...
%    saveimages, method)
%
% Variable Dictionary:
% --------------------
% imagename          input  The name of the image.
% imagefolder        input  The folder where the analysed images will be
%                           saved. This is the same as the input directory.
% vol                input  A 3D matrix version of the input TIFF image 
%                           stack/sequence. This is created by the
%                           PhantomAnalysisGUI using the loaded data.
% plane              input  The MRI imaging plane the image was acquired.
% hpixsize           input  The size of the voxel in X direction.
% vpixsize           input  The size of the voxel in Y direction.
% slicethickness     input  The size of the voxel in Z (slice) direction.
% centralslice       input  The central slice of the 3D matrix.
% ni                 input  The first slice to analyse.
% nf                 input  The last slice to analyse.
% thresholdf         input  Factor to multiply Canny detector's thresholds.
% sigma              input  Standard deviation of Canny detector's Gaussian
%                           blurring.
% saveimages         input  Equals 1 if user wants to save analysed images
%                           and 0 if not.
% method             input  The method for segmenting the phantom: 'Canny',
%                           'kmeans2' or 'kmeans3'.
% analysedvolmasks   output A 3D matrix of the analysed masks.
% volume             output The estimated volume (in mm^3).
%
% Last Modified: 02 February 2016
% Copyright (c) 2016, Xenios Milidonis

tic % stopwatch measuring elapsed time starts here

% Create new folder to save analysed images, if required.
if saveimages == 1
    mkdir(imagefolder,...
        sprintf('%s_phantomvolume', imagename(1:end - 4)));
    savefolder = [imagefolder,...
        sprintf('%s_phantomvolume', imagename(1:end - 4))];
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

% Get the position of the centroid and round it up. Using the un-scaled
% central slice is faster.
[centroidr, centroidc] = centroid(bwcentral);
centroidr = round(centroidr);
centroidc = round(centroidc);
  
% Call savevolume.m to save analysed images, if required.
if saveimages == 1
    savevolume(savefolder, centroidr, centroidc, central);
end

% Show in MATLAB's Command Window the headers for the measurements.
formatspec = '%9s %15s \n';
fprintf(formatspec, 'Slice', 'Volume');

% Create a column array to store the measurements.
volArray = zeros(nf - ni + 1, 1);

% Create 3 marker images to be used for reconstruction:
% a. one with a vertical white line through the centroid
% b. one with a white border of 2 pixels width around the image
% c. one with a white 5x5 area around the centroid ( for coronal images)
% For k-means check if images have a strong bias field by comparing the 
% mean of the intensity in the upper and lower half of the central slice.
% if strcmp(rotate, 'no')
    markervertical = zeros(size(central));
    markervertical(:, centroidc) = 1;
    markerborder = ones(size(central));
    markerborder(3:size(central, 1) - 2, 3:size(central, 2) - 2) = 0;  
    markercentroid = zeros(size(central));
    markercentroid(centroidr - 2:centroidr + 2,...
        centroidc - 2:centroidc + 2) = 1;

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
%     markervertical = zeros(size(rotated));
%     markervertical(:, centroidc) = 1;
%     markerborder = ones(size(rotated));
%     markerborder(3:size(rotated, 1) - 2, 3:size(rotated, 2) - 2) = 0; 
%     markercentroid = zeros(size(rotated));
%     markercentroid(centroidr - 2:centroidr + 2,...
%         centroidc - 2:centroidc + 2) = 1;
% 
%     meanAbove =...
%         mean(mean(rotated(centroidr - 50:centroidr, :)));
%     meanBelow =...
%         mean(mean(rotated(centroidr:centroidr + 50, :)));   
% end

% Main loop for measuring the volume of the phantom's central compartment.
for i = ni:nf
    rotated = vol(:, :, i);
    
%     % Rotate images by 90 degrees if required. For this case the centroid
%     % was re-calculated earlier.
%     if strcmp(rotate, 'yes')
%         rotated = imrotate(slice, -90);        
%     else
%         rotated = slice;
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
            % For Canny edge detection an iterated procedure is performed
            % by changing the size of Canny's Gaussian smoothing kernel
            % until the the segmented region fulfils a shape factor
            % criterion. The starting sigma is preserved for the next
            % slice.
            shapef = 0;
            sigmaloop = sigma;
            while shapef < 0.4  % or 0.35
                % Apply edge detection.
                [~, t] = edge(rotated, 'canny');
                [segmented, ~] = edge(rotated, 'canny', t * thresholdf,...
                    sigmaloop);
                
                % Remove some of the unwanted pixels only at the border of
                % the image. Internal pixels might help define the 
                % segmented areas.
                segmented = +segmented;
                boundarypixels = imreconstruct(markerborder, segmented);
                cleaned = logical(segmented - boundarypixels);

                % Fill the detected edges to create segmented regions. If
                % the image was acquired at the coronal/coronal90 plane do
                % a flood-fill at the centroid only.
                if strcmp(plane, 'Coronal') || strcmp(plane, 'Coronal90')
                    filled = fillandsmooth(cleaned, centroidr,...
                        centroidc, 'centroid');
                else
                    filled = fillandsmooth(cleaned, centroidr,...
                        centroidc, 'holes');
                end

                % Keep only the objects touching the vertical marker.
                filled = +filled;
                reconstructed = imreconstruct(markervertical, filled);

                % Keep only the largest object.
                largest = largestobject(reconstructed); 
                
                % Fill one last time to close any empty spots.
                largest = imfill(largest, 'holes'); 

                % Increase the size of the Gaussian smoothing kernel.
                sigmaloop = sigmaloop + 0.1;
                
                % Get the new shape factor. If < 0.37 repeat the analysis.
                shapef = shapefactor(largest, 'c');
                
                % Get the final analysed slice.
                analysedvolslice = +largest; 
            end
            
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

    % Post-processing steps for k-means.
    if strcmp(method, 'kmeans2') || strcmp(method, 'kmeans3') 
        % Remove some of the unwanted pixels at the border of the image and
        % internal pixels.
        segmented = +segmented;
        boundarypixels = imreconstruct(markerborder, segmented);
        cleaned = logical(segmented - boundarypixels);
        cleaned = bwareaopen(cleaned, ceil(size(cleaned, 1) / 7));

        % The regions are already filled. Keep only the objects touching
        % the vertical marker.
        cleaned = +cleaned;
        reconstructed = imreconstruct(markervertical, cleaned);

        % Keep only the largest object, or the one over the centroid.
        if strcmp(plane, 'Coronal') || strcmp(plane, 'Coronal90')    
            largest = imreconstruct(markercentroid, reconstructed);
        else
            largest = largestobject(reconstructed);
        end

        % Smooth the boundary of the segmented regions.
        filled = fillandsmooth(largest, centroidr, centroidc, 'holes');

        % Get the final analysed slice.
        analysedvolslice = +filled;
    end
    
    % Create a 3D matrix version of the analysed slices. This will be
    % passed as a function output to be shown in the GUI's 'Output Data'
    % figure.
    if i == ni                                     
        analysedvolmasks = analysedvolslice;
    else
        analysedvolmasks = cat(3, analysedvolmasks, analysedvolslice);
    end          
    
    % Count the number of pixels in each slice.
    if strcmp(method, 'Canny')
        slicesum = sum(sum(largest));
    else
        slicesum = sum(sum(filled));
    end
    
    % Fill array with counts.
    if isempty(slicesum) == 1
        volArray(i) = 0;
    else
        volArray(i) = slicesum;
    end
    
    % Call savevolumeloop.m to save analysed slices, if required.
    if saveimages == 1
        savevolumeloop(savefolder, rotated, segmented, cleaned,...
            filled, reconstructed, largest, i, centralslice, method)
    end 

    % Show in MATLAB's Command Window the measurements.
    formatspec = '%9d %15d\n';
    fprintf(formatspec, i, volArray(i));
end

% Calculate the total volume.
volumeinpixels = sum(volArray);
volume = volumeinpixels * hpixsize * vpixsize * slicethickness;

% Show in MATLAB's Command Window the total volume.
fprintf('%9s %15d\n', 'Total', volumeinpixels);
fprintf('%9s %15.3f\n', 'In mm^3', volume);

toc % stopwatch measuring elapsed time stops here

