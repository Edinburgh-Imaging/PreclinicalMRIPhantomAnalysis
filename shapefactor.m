function sf = shapefactor(bwimage, method)
% Companion function of the PhantomAnalysisGUI. It estimates the shape 
% factor of objects in a binary image using the formula of circularity 
% (4*pi*area / perimeter^2). There are 3 methods available for estimating
% the perimeter of objects:
% a. Count the boundary pixels.
% b. Sum the length of the sides of the pixels that lie on the boundary, 
%    given that a side has a length of 1. This is the best option for 
%    squarish objects.
% c. Sum the distance between each adjoining pair of pixels around the
%    boundary (default MATLAB method). This is the best option for rounded 
%    objects.
%
% >> sf = shapefactor(bwimage, method)
%
% Variable Dictionary:
% --------------------
% bwimage   input     The binary image.
% method    input     The method for calculating the perimeter: 'a', 'b' or
%                     'c'.
% sf        output    The estimated shape factor.
%
% Last Modified: 29 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Convert image to logical.
bwimage = logical(bwimage);

% Estimate the area of the object/s.
area = sum(sum(bwimage));

% Estimate the perimeter of the object/s.
switch method
    case 'a'  % METHOD a:
        perimeter = sum(sum(bwperim(bwimage)));
        
    case 'b'  % METHOD b:
        % Estimate the perimeter based on the sides of the pixels along the
        % boundary: get the complement of the image, with a border around
        % it in case the object is at the boundary.
        notimg = true(size(bwimage)+2);
        notimg(2:end-1, 2:end-1) = ~bwimage;

        % Find locations where a non-zero pixel is adjacent to a zero pixel
        % for each cardinal direction in turn.
        topedges =    bwimage & notimg(1:end-2, 2:end-1);
        leftedges =   bwimage & notimg(2:end-1, 1:end-2);
        bottomedges = bwimage & notimg(3:end,   2:end-1);
        rightedges =  bwimage & notimg(2:end-1, 3:end);

        % Sum each set of locations separately, then add to get total 
        % perimeter.
        perimeter = sum(topedges(:)) + sum(leftedges(:)) + ...
            + sum(bottomedges(:)) + sum(rightedges(:));
        
    case 'c'  % METHOD c:
        % Use regionprops to estimate the perimeter for each object and
        % list it in a structure.
        perimeterstruct = regionprops(bwimage, 'Perimeter');
        
        % Sum the perimeters, but first convert the structure to an array.
        perimeter = sum(cell2mat(struct2cell(perimeterstruct)));
end

% Estimate the shape factor of the object/s.
sf = 4 * pi * area / (perimeter ^ 2);

