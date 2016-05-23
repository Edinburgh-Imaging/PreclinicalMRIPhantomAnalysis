function bwattached = attachobjects(bwimage)
% Companion function of the PhantomAnalysisGUI. It attempts to attach
% multiple objects in the input binary image.
%
% >> bwattached = attachobjects(bwimage)
%
% Variable Dictionary:
% --------------------
% bwimage      input    The binary image.
% bwattached   output   The binary image with attached objects.
%
% Last Modified: 01 February 2016
% Copyright (c) 2016, Xenios Milidonis

% Count the number of white objects in the image.
[~, numofobjects] = bwlabel(bwimage);

% Morphologically open to remove small objects and then close to attach the
% remaining large objects.
if numofobjects > 1
    % Use a structural element based on the size of the image.
    se = strel('disk', round(size(bwimage, 1) / 30));
    bwimage = imopen(bwimage, se);
    bwattached = imclose(bwimage, se);
else
    bwattached = bwimage;
end

