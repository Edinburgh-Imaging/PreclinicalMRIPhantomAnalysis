function alignangle = phantomorientation(image, method)
% Companion function of the PhantomAnalysisGUI. It binarises the input
% image and then estimates an angle required to rotate and align the image
% with the horizontal x-axis. There are 2 methods available for estimating
% the angle of rotation:
% a. Fits an ellipse to the object and measures the angle between x-axis
%    and ellipse's major axis (default MATLAB method). This is the best
%    option for objects with perpendicular axes with differing lengths
%    (ellipse, rectangle etc).
% b. Identifies two points on the boundary of the object on each of the
%    left and right sides, the angle defined by each pair is measured and
%    the two angles are averaged. This is the best option for objects with
%    straight sides (square, rectangle etc).
%
% >> alignangle = phantomorientation(image, method)

% Variable Dictionary:
% --------------------
% image         input    The greyscale or binary image.
% method        input    The method for calculating the angle of rotation:
%                        'a' or 'b'.
% alignangle    output   The angle needed to align the image.
%
% Last Modified: 01 February 2016
% Copyright (c) 2016, Xenios Milidonis

% Binarise using k-means clustering (kmeanssegmentation.m). Two white
% clusters out of a total of 3 is the best for phantom images with no  
% uniform illumination. If the image is already a binary image, nothing
% changes.
bwimage = kmeanssegmentation(image, 3, 2);

% If multiple objects are found, try to attach them.
bwattached = attachobjects(bwimage);

% If there are still multiple objects, keep only the largest one.
bwlargest = largestobject(bwattached);

% Estimate the angle of rotation relative to the horizontal x-axis.
switch method
    case 'a'  % METHOD a:
        % Use regionprops to estimate the angle and list it in a structure.
        anglestruct = regionprops(bwlargest, 'Orientation');
        
        % The required angle to align the object is the inverse.
        alignangle = -anglestruct.Orientation;

    case 'b'  % METHOD b:
        % Get the position of the centroid and round it up.
        [centroidr, centroidc] = centroid(bwlargest);
        
        % Estimate the height of the object using a column vector.
        colvector = bwlargest(:, round(centroidc));
        height = sum(colvector);
        
        % Create 2 row vectors, one above and another below the centroid at
        % distances based on the height.
        rowup = round(centroidr - 0.3 * height);
        rowdo = round(centroidr + 0.3 * height);
        vectorup = bwlargest(rowup, :);
        vectordo = bwlargest(rowdo, :);
        
        % Find the column indices of the 4 points at the object boundary 
        % (2 left and 2 right).
        [~, colleftup]  = find(vectorup, 1, 'first');
        [~, colleftdo]  = find(vectordo, 1, 'first');        
        [~, colrightup] = find(vectorup, 1, 'last');
        [~, colrightdo] = find(vectordo, 1, 'last');  
        
        % Subtract down positions from up positions to make positions
        % relative to the origin (0, 0).
        rows = rowdo - rowup;
        colleft  = colleftdo  - colleftup;
        colright = colrightdo - colrightup;
        
        % The angle is the four-quadrant inverse tangent of the points.
        angleleft  = atan2d(rows, colleft);
        angleright = atan2d(rows, colright);
        
        % If the two angles are too different, ignore. Else, take their 
        % mean and subtract 90.
        if abs(angleleft - angleright) > 5
            alignangle = 0;
        else
            alignangle = mean([angleleft, angleright]) - 90;
        end
end

