function savedimensions(savefolder, ni, nf, centroidr, centroidc,...
    central, scaledseq)
% Companion function of the PhantomAnalysisGUI. It saves images pertaining
% to the analysis procedure for the estimation of the internal dimensions
% of the LEGO phantom before the main analysis loop.
%
% >> savedimensions(savefolder, ni, nf, centroidr, centroidc,...
%    central, scaledseq)
%
% Variable Dictionary:
% --------------------
% savefolder    input    The folder where the analysed images will be
%                        saved.
% ni            input    The first slice of the scaled matrix.
% nf            input    The last slice of the scaled matrix.
% centroidr     input    The row position of the centroid.
% centroidc     input    The column position of the centroid.
% central       input    The un-scaled central slice.
% scaledseq     input    The scaled 3D matrix.
%
% Last Modified: 27 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Get the central slice of the scaled matrix.
centralslicescaled = 1 + ceil((nf - ni) / 2);
centralscaled = scaledseq(:, :, centralslicescaled);

% STEP A. Save the central slice before analysis.
centralscaledpath = fullfile(savefolder, 'a_central_slice.png');
imwrite(centralscaled, centralscaledpath); 

% Call kmeanssegmentation.m to get the binarised central slice using 
% k-means clustering. Two white clusters out of a total of 3 is the best 
% option for phantom images with no uniform illumination.
bwcentral = kmeanssegmentation(central, 3, 2);

% STEP B. Save the thresholded central slice with a red cross indicating 
% the position of the centroid.
centralcentroid = insertMarker(bwcentral, [centroidc, centroidr],...
    'color', 'r', 'size', ceil(size(bwcentral, 1) / 30));
centralcentroidpath =...
    fullfile(savefolder, 'b_centroid_identification.png');
imwrite(centralcentroid, centralcentroidpath);

% STEP C. Save the thresholded central slice with a horizontal and a 
% vertical vector throught the centroid used for identification of the
% imaging plane.
[indexed, map] = gray2ind(bwcentral, 256);
rgbbwcentral = ind2rgb(indexed, map);
rgbbwcentral(centroidr, :, 1) = 0.7;
rgbbwcentral(:, centroidc, 1) = 0.7;
rgbbwcentralpath =...
    fullfile(savefolder, 'c_imaging_plane_identification.png');
imwrite(rgbbwcentral, rgbbwcentralpath);  
    
    
    

