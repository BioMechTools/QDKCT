function [oBoxSideL,oBox,boxMSR,boxes] = findBoundingBox2D( V,msr,minmax )
%Finds an extreme boundig box around a set of 2d point
%   Using the fact that at least one side of a minimum bounding box is
%   alligned with one side of the convex hull.
%
%Input:
%   V holds the points
%   msr is the the measure that is applied to the side length array
%       e.g. prod: the volume (default). sum: half the total side length 
%   minmax is either @min, or @max, which extreme box is selected
%   (default = @min)
%Output:
%   refL is the mean side length of the selected box
%   oBox holds the corner coordinates of the box
%       use plot(oBox([1:end,1],1),oBox([1:end,1],2))
%   boxMSR holds the values for the measure for all boxes
%   boxes holds all boxes
%
%       example use of these:
%       [~,minI] = min(boxMSR) ; [~,maxI] = max(boxMSR)  ;
%       cellfun(@(b) plot(b([1:end,1],1),b([1:end,1],2)),boxes([minI,maxI]),'uni',0) ;

if nargin < 2
    msr = @prod ;
end
if nargin < 3
    minmax = @min ;
end

[hullI] = convhull(V) ;
VwidRed = V(unique(hullI),:) ;
for i = numel(hullI)-1 : -1 : 1 ;
    ted = [hullI(i), hullI(i+1)] ;
    direction = V(ted(1),:)-V(ted(2),:) ;
    direction = direction / rssq(direction) ;
    %perpendicular vec in 2d found by crossing with an imaginary vec (3d)
    otherDir = cross([direction 0],[ 0 0 1]);
    otherDir = otherDir([1,2]) ;
    rotM = [direction;otherDir]  ;
    
    rotV= (rotM *VwidRed')' ;
    
    boxesCSL = [min(rotV);max(rotV)] ;
    
    if nargout > 1
        rotableBox = [boxesCSL([1 1 2 2 1],1),boxesCSL([1 2 2 1 1],2)] ;
        boxes{i} = (rotM\rotableBox')' ;
    end
    boxsideL = sort(abs(diff(boxesCSL))) ;
    boxMSR(i) = msr(boxsideL) ;
    boxSL(i,:) = boxsideL ;
end

[~,extremeI] = minmax(boxMSR) ;
% refL = mean(boxSL(extremeI,:)) ;
oBoxSideL =boxSL(extremeI,:) ;
if nargout >1
    oBox = boxes{extremeI} ;
end
end

 