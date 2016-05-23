function bwfilled = fillandsmooth(bwimage, centroidr, centroidc, where)
% Companion function of the PhantomAnalysisGUI. It uses morphological
% filling to define an area and smooth its boundary.
%
% >> bwfilled = fillandsmooth(bwimage, plane, centroidr, centroidc)
%
% Variable Dictionary:
% --------------------
% bwimage      input    The binary image.
% centroidr    input    The row position of the centroid.
% centroidc    input    The column position of the centroid.
% where        input    Where to fill: 'holes', 'centroid'.
% bwsmoothed   output   The smoothed binary image.
%
% Last Modified: 31 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Create two structural elements.
sev = strel('line', 5, 90);  % vertical
seh = strel('line', 3, 0);   % horizontal

% Dilate the image to connect adjacent pixels.
bwdilated = imdilate(bwimage, [sev, seh]);

% Fill the image at the centroid or at all holes.
if strcmp(where, 'centroid')
    bwfilled = imfill(bwdilated, [centroidr, centroidc]); 
else
    bwfilled = imfill(bwdilated, 'holes'); 
end

% Erode the image to bring the objects back to their correct size.
bwfilled = imerode(bwfilled, [sev, seh]);

