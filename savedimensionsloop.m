function savedimensionsloop(savefolder, ni, nf, centroidr, centroidc,...
    scalefactor, rotated, edged, cleaned, shiftUD, shiftUD2,...
    shiftLR, shiftLR2, i)
% Companion function of the PhantomAnalysisGUI. It saves images pertaining
% to the analysis procedure for the estimation of the internal dimensions
% of the LEGO phantom inside the main analysis loop, such as a masks for 
% each slice.
%
% >> savedimensionsloop(savefolder, ni, nf, centroidr, centroidc,...
%    scalefactor, rotated, edged, cleaned, shiftUD, shiftUD2,...
%    shiftLR, shiftLR2, i)
%
% Variable Dictionary:
% --------------------
% savefolder    input    The folder where the analysed images will be
%                        saved.
% ni            input    The first slice of the scaled matrix.
% nf            input    The last slice of the scaled matrix.
% centroidr     input    The row position of the centroid.
% centroidc     input    The column position of the centroid.
% scalefactor   input    The factor with which the matrix was scaled.
% rotated       input    The rotated central slice.
% edged         input    The boundary image of the central slice.
% cleaned       input    The cleaned version of 'edged'.
% shiftUD, shiftUD2, shiftLR, shiftLR2
%               input    Various distances from the centroid in both
%                        directions used to create vectors.
% i             input    The index of the slice. 
%
% Last Modified: 27 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Get the index of the central slice of the scaled matrix.
centralslicescaled = 1 + ceil((nf - ni) / 2);

% Further analyses on the central slice only.
if i == centralslicescaled
    % STEP D. Save the rotated central slice.
    rotatedpath = fullfile(savefolder, 'd_rotation.png');
    imwrite(rotated, rotatedpath);
    
    % STEP E. Save the boundary image of the central slice.
    edgedpath = fullfile(savefolder, 'e_boundary_detection.png');
    imwrite(edged, edgedpath);
    
    % STEP F. Save the cleaned version of 'edged'.
    cleanedpath = fullfile(savefolder, 'f_cleaning.png');
    imwrite(cleaned, cleanedpath); 
    
    % STEP G. Save the cleaned slice with a horizontal vector throught the
    % centroid used for providing a first rought estimate of the horizontal
    % phantom dimension.
    [indexed, map] = gray2ind(cleaned, 256);
    rgbcleaned = ind2rgb(indexed, map);
    rgbcleaned(round(centroidr * scalefactor), :, 1) = 0.7;
    rgbcleanedpath = fullfile(savefolder,...
        'g_rough_horizontal_dimension_estimation.png');
    imwrite(rgbcleaned, rgbcleanedpath);  
    
    % STEP H. Save the cleaned slice with all horizontal and vertical
    % vectored used for measuring the dimensions of the phantom.
    [indexedvectors, map] = gray2ind(cleaned, 256);
    rgbvectors = ind2rgb(indexedvectors, map);
    rgbvectors(round(centroidr * scalefactor) + shiftUD, :, 1) = 0.7;
    rgbvectors(round(centroidr * scalefactor) - shiftUD, :, 1) = 0.7;
    rgbvectors(:, round(centroidc * scalefactor) + shiftLR, 1) = 0.7;
    rgbvectors(:, round(centroidc * scalefactor) - shiftLR, 1) = 0.7;
    rgbvectors(round(centroidr * scalefactor) + shiftUD2, :, 1) = 0.7;
    rgbvectors(round(centroidr * scalefactor) - shiftUD2, :, 1) = 0.7;
    rgbvectors(:, round(centroidc * scalefactor) + shiftLR2, 1) = 0.7;
    rgbvectors(:, round(centroidc * scalefactor) - shiftLR2, 1) = 0.7;
    vectorspath = fullfile(savefolder, 'h_dimensions_estimation.png');
    imwrite(rgbvectors, vectorspath);
end

% Save a mask for each slice.
icleanedpath = fullfile(savefolder, sprintf('%d.png', i + ni - 1));
imwrite(cleaned, icleanedpath);  

% Save each rotated slice (before other analyses) with edges overlaid.
[indexed, map] = gray2ind(rotated, 256);
rgboverlaid = ind2rgb(indexed, map);
[rows, columns] = find(cleaned);
for j = 1:length(rows)
    rgboverlaid(rows(j), columns(j), 1) = 0.7; % 'burns' the pixels
end 
rgboverlaidpath = fullfile(savefolder,...
    sprintf('%d_overlaid.png', i + ni - 1));
imwrite(rgboverlaid, rgboverlaidpath);

