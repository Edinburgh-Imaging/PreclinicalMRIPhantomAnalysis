function [horuArray, hordArray, verlArray, verrArray, horuArray2,...
    hordArray2, verlArray2, verrArray2, locationsArray, shiftUD,...
    shiftUD2, shiftLR, shiftLR2] =...
    vectorsdimensions(cleaned, centroidr, centroidc, scalefactor,...
    horuArray, hordArray, verlArray, verrArray, horuArray2, hordArray2,...
    verlArray2, verrArray2, locationsArray, plane, hpixsize, vpixsize,...
    horpixels, i, method)
% Companion function of the PhantomAnalysisGUI. It creates horizontal and
% vertical vectors at positions relative to the centroid of the phantom
% that are used for estimating the internal dimensions of the phantom.
%
% >> [horuArray, hordArray, verlArray, verrArray, horuArray2,...
%    hordArray2, verlArray2, verrArray2, locationsArray, shiftUD,...
%    shiftUD2, shiftLR, shiftLR2] =...
%    vectorsdimensions(cleaned, centroidr, centroidc, scalefactor,...
%    horuArray, hordArray, verlArray, verrArray, horuArray2, hordArray2,...
%    verlArray2, verrArray2, locationsArray, plane, hpixsize, vpixsize,... 
%    horpixels, i, method)
%
% Variable Dictionary:
% --------------------
% cleaned        input   The cleaned version of 'edged'.
% centroidr      input   The row position of the centroid.
% centroidc      input   The column position of the centroid.
% scalefactor    input   The factor with which the matrix was scaled.
% horuArray, hordArray, verlArray, verrArray, horuArray2, hordArray2,...
% verlArray2, verrArray2
%                input   Empty column arrays for the measurements.
% locationsArray input   A matrix containing the locations where vectors 
%                        and edges overlap.
% plane          input   The MRI imaging plane the image was acquired.
% hpixsize       input   The size of the voxel in X direction.
% vpixsize       input   The size of the voxel in Y direction.
% horpixels      input   The rough horizontal dimension of the phantom.
% i              input   The index of the slice.
% method         input   'Canny', 'kmeans2' or 'kmeans3'. Required 
%                        identifier as in k-means clustering 1 pixel 
%                        distance must be included in the calculation to 
%                        take into account both pixels of the identified 
%                        boundary.
% horuArray, hordArray, verlArray, verrArray, horuArray2, hordArray2,...
% verlArray2, verrArray2
%                output  Filled column arrays with the measurements.
% shiftUD, shiftUD2, shiftLR, shiftLR2
%                output  Various distances from the centroid in both
%                        directions used to create the vectors.
%
% Last Modified: 31 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Based on the horizontal dimension, create 4 horizontal vectors shifted
% up and down from centroid.
shiftUD  = ceil(horpixels / 13); 
shiftUD2 = ceil(horpixels / 5); 
horvectorupos  = round(centroidr * scalefactor) + shiftUD;
horvectordpos  = round(centroidr * scalefactor) - shiftUD;
horvectoru2pos = round(centroidr * scalefactor) + shiftUD2;
horvectord2pos = round(centroidr * scalefactor) - shiftUD2;
horvectoru  = cleaned(horvectorupos,  :);
horvectord  = cleaned(horvectordpos,  :); 
horvectoru2 = cleaned(horvectoru2pos, :);
horvectord2 = cleaned(horvectord2pos, :); 

% Find the positions where the vectors overlap the phantom boundary and
% estimate the distance between opposing pairs.
[~, horleftcolu]   = find(horvectoru,  1, 'first');
[~, horrightcolu]  = find(horvectoru,  1, 'last');
[~, horleftcold]   = find(horvectord,  1, 'first');
[~, horrightcold]  = find(horvectord,  1, 'last');                  
[~, horleftcolu2]  = find(horvectoru2, 1, 'first');
[~, horrightcolu2] = find(horvectoru2, 1, 'last');
[~, horleftcold2]  = find(horvectord2, 1, 'first');
[~, horrightcold2] = find(horvectord2, 1, 'last');
if strcmp(method, 'Canny')
    horpixelsu  = horrightcolu  - horleftcolu;
    horpixelsd  = horrightcold  - horleftcold;
    horpixelsu2 = horrightcolu2 - horleftcolu2;
    horpixelsd2 = horrightcold2 - horleftcold2;
else
    horpixelsu  = horrightcolu  - horleftcolu  + 1;
    horpixelsd  = horrightcold  - horleftcold  + 1;
    horpixelsu2 = horrightcolu2 - horleftcolu2 + 1;
    horpixelsd2 = horrightcold2 - horleftcold2 + 1;  
end
hordistanceu  = horpixelsu  * hpixsize / scalefactor;
hordistanced  = horpixelsd  * hpixsize / scalefactor; 
hordistanceu2 = horpixelsu2 * hpixsize / scalefactor;
hordistanced2 = horpixelsd2 * hpixsize / scalefactor;
 

% Based on the horizontal dimension, create 4 vertical vectors shifted
% left and right from centroid. The position of the vectors for axial/
% axial90/sagittal/sagittal90 planes will be different.
if strcmp(plane, 'Coronal') || strcmp(plane, 'Coronal90')
    shiftLR  = ceil(horpixels / 13);
    shiftLR2 = ceil(horpixels / 5);
else
    shiftLR  = ceil(horpixels * 3 / 9);
    shiftLR2 = ceil(horpixels * 4 / 9);
end
vervectorlpos  = round(centroidc * scalefactor) - shiftLR;
vervectorrpos  = round(centroidc * scalefactor) + shiftLR;
vervectorl2pos = round(centroidc * scalefactor) - shiftLR2;
vervectorr2pos = round(centroidc * scalefactor) + shiftLR2;
vervectorl  = cleaned(:, vervectorlpos);
vervectorr  = cleaned(:, vervectorrpos);
vervectorl2 = cleaned(:, vervectorl2pos);
vervectorr2 = cleaned(:, vervectorr2pos);

% Find the positions where the vectors overlap the phantom boundary and
% estimate the distance between opposing pairs.
[verupperrowl, ~]  = find(vervectorl,  1, 'first');
[verlowerrowl, ~]  = find(vervectorl,  1, 'last');
[verupperrowr, ~]  = find(vervectorr,  1, 'first');
[verlowerrowr, ~]  = find(vervectorr,  1, 'last');
[verupperrowl2, ~] = find(vervectorl2, 1, 'first');
[verlowerrowl2, ~] = find(vervectorl2, 1, 'last');
[verupperrowr2, ~] = find(vervectorr2, 1, 'first');
[verlowerrowr2, ~] = find(vervectorr2, 1, 'last');
if strcmp(method, 'Canny')
    verpixelsl  = verlowerrowl  - verupperrowl;
    verpixelsr  = verlowerrowr  - verupperrowr;
    verpixelsl2 = verlowerrowl2 - verupperrowl2;
    verpixelsr2 = verlowerrowr2 - verupperrowr2;
else
    verpixelsl  = verlowerrowl  - verupperrowl  + 1;
    verpixelsr  = verlowerrowr  - verupperrowr  + 1;
    verpixelsl2 = verlowerrowl2 - verupperrowl2 + 1;
    verpixelsr2 = verlowerrowr2 - verupperrowr2 + 1;
end
verdistancel  = verpixelsl  * vpixsize / scalefactor;
verdistancer  = verpixelsr  * vpixsize / scalefactor;        
verdistancel2 = verpixelsl2 * vpixsize / scalefactor;
verdistancer2 = verpixelsr2 * vpixsize / scalefactor; 

% Fill the empty column arrays with the measurements. Also, save the 
% positions where the vectors overlap the boundary of the phantom to be
% used for visualisation in the GUI's 'Output Data' figure. These should be
% projected back to the original size of the image.
if isempty(hordistanceu) == 1
    horuArray(i) = 0;
else
    horuArray(i) = hordistanceu;
    locationsArray(1, :, i) = [round(horvectorupos / scalefactor),...
        round(horleftcolu / scalefactor)];
    locationsArray(2, :, i) = [round(horvectorupos / scalefactor),...
        round(horrightcolu / scalefactor)];
end
if isempty(hordistanced) == 1
    hordArray(i) = 0;
else
    hordArray(i) = hordistanced;
    locationsArray(3, :, i) = [round(horvectordpos / scalefactor),...
        round(horleftcold / scalefactor)];
    locationsArray(4, :, i) = [round(horvectordpos / scalefactor),...
        round(horrightcold / scalefactor)];        
end
if isempty(hordistanceu2) == 1
    horuArray2(i) = 0;
else
    horuArray2(i) = hordistanceu2;
    locationsArray(5, :, i) = [round(horvectoru2pos / scalefactor),...
        round(horleftcolu2 / scalefactor)];
    locationsArray(6, :, i) = [round(horvectoru2pos / scalefactor),...
        round(horrightcolu2 / scalefactor)];
end
if isempty(hordistanced2) == 1
    hordArray2(i) = 0;
else
    hordArray2(i) = hordistanced2;
    locationsArray(7, :, i) = [round(horvectord2pos / scalefactor),...
        round(horleftcold2 / scalefactor)];
    locationsArray(8, :, i) = [round(horvectord2pos / scalefactor),...
        round(horrightcold2 / scalefactor)]; 
end
if isempty(verdistancel) == 1
    verlArray(i) = 0;
else
    verlArray(i) = verdistancel;
    locationsArray(9, :, i) = [round(verupperrowl / scalefactor),...
        round(vervectorlpos / scalefactor)];
    locationsArray(10, :, i) = [round(verlowerrowl / scalefactor),...
        round(vervectorlpos / scalefactor)];
end
if isempty(verdistancer) == 1
    verrArray(i) = 0;
else
    verrArray(i) = verdistancer;
    locationsArray(11, :, i) = [round(verupperrowr / scalefactor),...
        round(vervectorrpos / scalefactor)];
    locationsArray(12, :, i) = [round(verlowerrowr / scalefactor),...
        round(vervectorrpos / scalefactor)];
end
if isempty(verdistancel2) == 1
    verlArray2(i) = 0;
else
    verlArray2(i) = verdistancel2;
    locationsArray(13, :, i) = [round(verupperrowl2 / scalefactor),...
        round(vervectorl2pos / scalefactor)];
    locationsArray(14, :, i) = [round(verlowerrowl2 / scalefactor),...
        round(vervectorl2pos / scalefactor)];
end
if isempty(verdistancer2) == 1
    verrArray2(i) = 0;
else
    verrArray2(i) = verdistancer2;
    locationsArray(15, :, i) = [round(verupperrowr2 / scalefactor),...
        round(vervectorr2pos / scalefactor)];
    locationsArray(16, :, i) = [round(verlowerrowr2 / scalefactor),...
        round(vervectorr2pos / scalefactor)];
end       
    
    