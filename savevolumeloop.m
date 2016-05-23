function savevolumeloop(savefolder, rotated, segmented, cleaned,...
    filled, reconstructed, largest, i, centralslice, method)
% Companion function of the PhantomAnalysisGUI. It saves images pertaining
% to the analysis procedure for the estimation of the volume of the central
% frustum-shaped compartment of the LEGO phantom inside the main analysis 
% loop, such as a masks for each slice.
%
% >> savevolumeloop(savefolder, rotated, segmented, cleaned,...
%    filled, bwreconstructed, bwlargest, i, centralslice)
%
% Variable Dictionary:
% --------------------
% savefolder    input    The folder where the analysed images will be
%                        saved.
% rotated       input    The rotated central slice.
% segmented     input    The segmented image of the central slice.
% cleaned       input    The cleaned version of 'edged'.
% filled        input    The filled image.
% reconstructed input    The reconstructed image.
% largest       input    The image with the largest object.
% i             input    The index of the slice. 
% centralslice  input    The central slice of the 3D matrix.
% method        input    The method for segmenting the phantom: 'Canny',
%                        'kmeans2' or 'kmeans3'.
%
% Last Modified: 27 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Further analyses on the central slice only. The order is different for
% Canny and k-means.
if i == centralslice
    % STEP D. Save the rotated central slice.
    rotatedpath = fullfile(savefolder, 'd_rotation.png');
    imwrite(rotated, rotatedpath);
    
    % STEP E. Save the boundary image of the central slice.
    segmentedpath = fullfile(savefolder, 'e_segmentation.png');
    imwrite(segmented, segmentedpath);
    
    % STEP F. Save the cleaned version of 'bwimage'.
    cleanedpath = fullfile(savefolder, 'f_cleaning.png');
    imwrite(cleaned, cleanedpath); 

    if strcmp(method, 'Canny')
        % STEP G. Save the filled image.
        filledpath = fullfile(savefolder, 'g_filling.png');
        imwrite(filled, filledpath);  

        % STEP H. Save the reconstructed image touching the vertical mask.
        reconstructedpath = fullfile(savefolder, 'h_reconstructing.png');
        imwrite(reconstructed, reconstructedpath);  

        % STEP I. Save the image with the largest object only.
        largestpath = fullfile(savefolder, 'i_largest.png');
        imwrite(largest, largestpath);          
    else
        % STEP G. Save the reconstructed image touching the vertical mask.
        reconstructedpath = fullfile(savefolder, 'g_reconstructing.png');
        imwrite(reconstructed, reconstructedpath);  

        % STEP H. Save the image with the largest object only.
        largestpath = fullfile(savefolder, 'h_largest.png');
        imwrite(largest, largestpath);  
    
        % STEP I. Save the filled image.
        filledpath = fullfile(savefolder, 'i_filling.png');
        imwrite(filled, filledpath);    
    end
end

% Save a mask for each slice.
if strcmp(method, 'Canny')
    ilargestpath = fullfile(savefolder, sprintf('%d.png', i));
    imwrite(largest, ilargestpath);
else
    ifilledpath = fullfile(savefolder, sprintf('%d.png', i));
    imwrite(filled, ifilledpath);
end

% Save each rotated slice (before other analyses) with edges overlaid.
[indexed, map] = gray2ind(rotated, 256);
rgboverlaid = ind2rgb(indexed, map); 
if strcmp(method, 'Canny')
    [rows, columns] = find(largest); 
else
    [rows, columns] = find(filled);
end
for j = 1:length(rows)
    rgboverlaid(rows(j), columns(j), 1) = 0.7; % 'burns' the pixels
end 
rgboverlaidpath = fullfile(savefolder, sprintf('%d_overlaid.png', i));
imwrite(rgboverlaid, rgboverlaidpath);

