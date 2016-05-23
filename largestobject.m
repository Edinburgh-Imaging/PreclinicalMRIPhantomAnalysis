function bwlargest = largestobject(bwimage)
% Companion function of the PhantomAnalysisGUI. It identifies and keeps
% only the largest object in a binary image.
%
% >> bwlargest = largestobject(bwimage)
%
% Variable Dictionary:
% --------------------
% bwimage      input    The binary image with multiple objects.
% bwlargest    output   The binary image with only the largest object.
%
% Last Modified: 27 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Convert image to binary.
bwimage = logical(bwimage);

% Create an empty image with the same size.
bwlargest = zeros(size(bwimage));

% Get the area and index of each object in the image.
areastruct = regionprops(bwimage, 'Area', 'PixelIdxList');

% Identify the object with the largest area.
[~, index] = max([areastruct.Area]);

% Copy the largest object into the empty image.
bwlargest(areastruct(index).PixelIdxList) = 1;