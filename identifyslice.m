function [firstslice, centralslice, lastslice] = identifyslice(vol)
% Companion function of the PhantomAnalysisGUI. It identifies the first, 
% central and last slices through the phantom and phantom's central 
% frustum-shaped compartment (in z direction).
%
% >> [firstslice, centralslice, lastslice] = identifyslice(vol)
%
% Variable Dictionary:
% --------------------
% vol           input    The greyscale 3D matrix.
% firstslice    output   The estimated number of the first slice.
% centralslice  output   The estimated number of the central slice.
% lastslice     output   The estimated number of the last slice.
%
% Last Modified: 29 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Normalise image intensities from 0 to 1.
vol = vol - min(vol(:));
vol = vol / max(vol(:));

% Smooth each slice and create a new volume.
for i = 1:size(vol, 3)
    h = fspecial('gaussian', [3, 3], 1);
    blurred = imfilter(vol(:, :, i), h);

    if i == 1                                       
        smoothedvol = blurred;
    else
        smoothedvol = cat(3, smoothedvol, blurred);
    end
end

% Call kmeanssegmentation.m to get the binarised volume slice using 
% k-means clustering. One white cluster and one black cluster is the best 
% option for data with a large background noise distribution. 
% This is almost identical to binarisation with Otsu's threshold.
bwvol = kmeanssegmentation(smoothedvol, 2, 1);

% FIND THE LARGEST CONNECTED COMPONENT IN 3D SPACE.
% Get all components (6 is the least dense connectivity).
conncomp = bwconncomp(bwvol, 6);   

% Identify the largest component using cellfun.
[~, maxcell] = max(cellfun(@numel, conncomp.PixelIdxList));  

% Zero the mask and assign to it the largest component.
bwvol = zeros(size(bwvol));
bwvol(conncomp.PixelIdxList{1, maxcell}) = 1;

% Remove smaller regions in each slice in an attempt to exclude partial
% volume effects at the edges, based on the total volume, and identify the
% slices that include some of the object.
totalpixels = sum(bwvol(:));
objectindex = zeros(size(vol, 3), 1);
for i = 1:size(vol, 3)
    cleaned = bwareaopen(bwvol(:, :, i), round(totalpixels/30));

    % Fill empty array with 1 where the object is present.
    if nnz(cleaned) > 0
        objectindex(i) = 1;
    end                                    
end
 
% Estimate the central slice.
first = find(objectindex, 1, 'first');
last = find(objectindex, 1, 'last');
centralslice = ceil((last + first)/2);

% Estimate the first and last slices using the slice range of the central 
% phantom compartment 'centralslices'.
totalslices = last - first + 1;
centralslices = totalslices * 5 / 12.8;
firstslice = centralslice - floor(centralslices/2);
lastslice = centralslice + floor(centralslices/2);

% If this fails, just estimate the central slice of the whole 3D matrix.
if isempty(centralslice) == 1
    centralslice = ceil(size(vol, 3)/2);
    firstslice = centralslice - 2;
    lastslice = last + 2;
end

