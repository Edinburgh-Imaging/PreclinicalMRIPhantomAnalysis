function [centroidr, centroidc] = centroid(bwimage)
% Companion function of the PhantomAnalysisGUI. It estimates the position
% of the centroid of a binary image. The centroid is estimated by finding 
% the mean of the positions of all white pixels in the image in the two 
% directions. This provides identical estimates to MATLAB's 
% regionprops(bwimage, 'Centroid') method.
%
% >> [centroidr, centroidc] = centroid(bwimage)
%
% Variable Dictionary:
% --------------------
% bwimage      input    The binary image.
% centroidr    output   The estimated row of the centroid.
% centroidc    output   The estimated column of the centroid.
%
% Last Modified: 01 February 2016
% Copyright (c) 2016, Xenios Milidonis

% Find the positions of all white pixels.
[rowsofwhite, colsofwhite] = find(bwimage);

% The centroid is the mean of horizontal and vertical positions.
centroidr = mean(rowsofwhite);
centroidc = mean(colsofwhite);

