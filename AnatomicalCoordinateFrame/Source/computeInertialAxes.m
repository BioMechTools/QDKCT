function  [CM,IA,SV,V,F,VcsL] = computeInertialAxes(F,V,gridSize,percRange,Axes)
%% COMPUTEINERTIALAXES(F,V,gridSize[,percRange,Axes]) Determine the inertial axes of a surface mesh by filling the volume uniformly with point-masses on a cubic grid.
%[CM,IA,SV,V,F,VcsL] = computeInertialAxes(F,V,gridSize,percRange,Axes)
%
%---Input
%     F            -    Faces (N by 3 integer
%     V            -    Vertices (N by 3)
%     gridSize     -    size of the point-mass grid (use same units as V)
%     [percRange]  -    Optionally only a section (in percentage range) can
%            be evaluated
%     [Axes]       -    The percentage range is taken along the direction
%     of the first ROW of the reference frame given in Axes
%   
%---Output
%     CM           -     Center of volume mass (1 by 3 double)
%     IA           -     Inertial axes (3x3, each column is a direction)
%     SV           -     Corresponding moments of inertia
%       When percRange and Axes are supplied: 
%     V            -     Vertices of the section
%     F            -     Faces of the section
%     VcsL         -     Vertices of the section in Axes coordinate system
%
%  see also:  cutMeshByPercent 

if nargin > 3
    rotM = eye(3) * Axes' ;
    V = (rotM * V')' ;
    [F, V] = cutMeshByPercent(F,V,percRange,1,true);
end
if nargin <3 || isempty(gridSize)
    gridSize = 0.5;
end

[OV, OVpos] = surface2volume_nested(F,V,[],0,gridSize);
[CM,IA,SV] = InertialAxesVolume_nested(OV,OVpos);

if nargin > 3
    IA = (rotM \ IA) ;
    CM = (rotM \ CM')' ;
    if nargout > 3
        if nargout > 5
            VcsL = V;
        end
        V = (rotM \ V')' ;
    end
end

[SV,sortI] = sort(SV,'ascend') ;
IA = IA(:,sortI) ;
end



function [outputVolume, outputGrid] = surface2volume_nested(Faces,Vertices,inputGrid,verboseOutput,gridSize)
%SURFACE2VOLUME convert a surface volume to a solid volume
%   [outputVolume,outputGrid] = surface2volume_nested(Faces,Vertices) creates a 
%   volume block (logical) in which every voxel which is inside the given 
%   surface is set to 1. The output block has the dimension [M,N,P] where 
%       M = size(min(Xpos) - 2*GRIDSIZE : GRIDSIZE : max(Xpos) + 2*GRIDSIZE)
%       N = size(min(Ypos) - 2*GRIDSIZE : GRIDSIZE : max(Ypos) + 2*GRIDSIZE)
%       P = size(min(Zpos) - 2*GRIDSIZE : GRIDSIZE : max(Zpos) + 2*GRIDSIZE)
%   Xpos, Ypos, Zpos are the Coordinates of the given vertices. The default 
%   GRIDSIZE is 0.5.
%
%   [outputVolume,outputGrid] = surface2volume_nested(Faces,Vertices,{X,Y,Z}) uses 
%   the grid coordinates given by X,Y,Z instead of creating an own grid. 
%   The grid coordinates have to be uniform distributed as given by the 
%   output of meshgrid using the same distance in all three directions 
%       [X,Y,Z] = meshgrid(1:GR:M,1:GR:N,1:GR:P);
%   The script tests slightly if the given grid is equidistant. The
%   definition of an input grid allows also to cut the surface, but this
%   should be used with care.
%   
%   X Y and Z have to be assigned as a cell array using curly braces.
%   
%   [outputVolume,outputGrid] = surface2volume_nested(Faces,Vertices,{X,Y,Z},1) 
%   writes information about the progress to the standard output.
%
%   How it works:
%       First the surface will be rasterized on the grid. Therefore it 
%       calculates the position of points which lie in the surface in a
%       finer resolution as defined by the inputgrid. These points were
%       then tranfered to the point it the inputgrid by using a simple
%       indexing technique. One could also use dsearchn, but this takes to
%       much computational time, however, it can avoid the need of
%       an equidistant grid. After rasterizing the patches the
%       background is fill using imfill. The start point is set to the
%       lower left corner. Afterwards the datablock will be inverted. The
%       script tests if the datablock is fully filled  and tries to repeat 
%       the task slice by slice.
%
%   Suggestions: Use a very fine surface. Sometimes it fails if the patches
%   of the surfaces are to large and the volume isn't closed after
%   rasterization.
%
%   Example:
%       load mri;
%       D = squeeze(D);
%       D = padarray(D,[5 5 5],'both');
%       
%       % Create an isosurface
%       Ds = smooth3(D);
%       surface = isosurface(Ds,5);
%       
%       % Display the surface
%       figure;
%       subplot(1,2,1);
%       hiso = patch('Vertices',Vertices,...
%                    'Faces',Faces,...
%                    'FaceColor',[1,.75,.65],...
%                    'EdgeColor','none');
%
%       view(45,30) 
%       axis tight 
%       daspect([1,1,.4])
%       lightangle(45,30); 
%       set(gcf,'Renderer','zbuffer'); lighting phong
%       isonormals(Ds,hiso)
%       set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)
%
%       % Reconstruct the volume and display it as montage
%       outputVolume = surface2volume_nested(Faces,Vertices,[],1);
%       nDims = size(outputVolume);
%       subplot(1,2,2);
%       montage(reshape(outputVolume,nDims(1),nDims(2),1,nDims(3)),[0 1]);

DEFAULT_GRIDSIZE = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%
% Check input function %
%%%%%%%%%%%%%%%%%%%%%%%%

% Report command window
if (~exist('verboseOutput','var')) 
    verboseOutput = 0;
end
outputFID = 1;
if verboseOutput
    fprintf(outputFID,'surface2volume_nested \n') ;
end

if verboseOutput
    tic; fprintf(outputFID,'Initializing ... ');
end
% Check surface consist of triangles
nFaces = size(Faces,1); 
if (size(Faces,2) ~= 3)
   error('Matlab:surface2volume_nested','Input faces must be triangles.'); 
end

if (~(exist('gridSize','var')))
    gridSize = DEFAULT_GRIDSIZE;
end

if (~exist('inputGrid','var'))
    inputGrid = [];
end

% inputGrid
if (isempty(inputGrid))
    minVPos = min(Vertices,[],1);
    maxVPos = max(Vertices,[],1);
    [Xgrid,Ygrid,Zgrid] = ndgrid( minVPos(1) - 2*gridSize : gridSize : maxVPos(1) + 2*gridSize, ...
                        minVPos(2) - 2*gridSize : gridSize : maxVPos(2) + 2*gridSize, ...
                        minVPos(3) - 2*gridSize : gridSize : maxVPos(3) + 2*gridSize);
    inputGrid = {Xgrid Ygrid Zgrid};
else
    % test for gridsize of inputgrid
    grX = diff(inputGrid{1,1}([1,2]));
    grY = diff(inputGrid{1,2}([1,2]));
    grZ = diff(inputGrid{1,3}([1,2]));
    if (~(grX || grY || grZ))
       error('Matlab:surface2volume_nested','The input grid is not equidistant.');
    else
        sizeV = [grX grY grZ];
        pV = find(sizeV);
        gridSize = sizeV(pV(1));
    end
end

startGridPosX = min(inputGrid{1,1}(:)); % start grid position
startGridPosY = min(inputGrid{1,2}(:));
startGridPosZ = min(inputGrid{1,3}(:));

endGridPosX = max(inputGrid{1,1}(:)); % end grid position
endGridPosY = max(inputGrid{1,2}(:));
endGridPosZ = max(inputGrid{1,3}(:));

if (~iscell(inputGrid))
    error('Matlab:surface2volume_nested','Input grid must be a cell array.');
end

outputGrid = inputGrid; % output when no inputGrid was defined

initialVSize = size(inputGrid{1,1});
outputVolume = zeros(initialVSize,'uint8'); % Dummy outputVolume

if verboseOutput
    fprintf(outputFID,'done in %g sec\n',toc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% raster faces to points %
%%%%%%%%%%%%%%%%%%%%%%%%%%
if verboseOutput
    tic; fprintf(outputFID,'Rasterize points in patches to grid points ... ');
end

for iFace = 1:nFaces
    pointMatrix = Vertices(Faces(iFace,:),:);
    % test if triangle is in the search volume 
    % (this cost around 10% performance, if were is nothing to omit, 
    % however it can reduce the amount of cycles, if the search volume 
    % is smaller than the surface object)
    inX = nnz(((pointMatrix(:,1) >= startGridPosX) & (pointMatrix(:,1) <= endGridPosX)));
    inY = nnz(((pointMatrix(:,2) >= startGridPosY) & (pointMatrix(:,2) <= endGridPosY)));
    inZ = nnz(((pointMatrix(:,3) >= startGridPosZ) & (pointMatrix(:,3) <= endGridPosZ)));
    if (inX  && inY && inZ)
        tPoints = pointsontriangle_nested(pointMatrix, gridSize/2);
        % round positions on grid indices
        tPoints(:,1) = round((tPoints(:,1) - startGridPosX + gridSize)./gridSize);
        tPoints(:,2) = round((tPoints(:,2) - startGridPosY + gridSize)./gridSize);
        tPoints(:,3) = round((tPoints(:,3) - startGridPosZ + gridSize)./gridSize);
        % set the indices belonging to the points in the outputvolume to 1
        nTestPoints = size(tPoints,1);
        for iTestPoint = 1:nTestPoints
            outputVolume(tPoints(iTestPoint,1),tPoints(iTestPoint,2),tPoints(iTestPoint,3)) = 1;
        end 
    end
end
% crop outputvolume if some point were found outside
if (size(outputVolume) ~= initialVSize);
    outputVolume = outputVolume(1:initialVSize(1),1:initialVSize(2),1:initialVSize(2));
end


if verboseOutput
    fprintf(outputFID,'done in %g sec\n',toc);
end
%%%%%%%%%%%%%%%
% fill volume %
%%%%%%%%%%%%%%%
% fill the background of the volume using imfill, assuming that the lower
% left corner (beginning of the block is located in background)
if verboseOutput
    tic; fprintf(outputFID,'Fill volume ... ');
end

volumeShape = outputVolume;

outputVolume = ~(imfill(logical(outputVolume),[1 1 1],6));
if (all(outputVolume) |  ~any(outputVolume))
    warning('Matlab:surface2volume_nested','The output is full flooded, assume that the surface was not close. Retrying slicewise.');
    for i=1:initialVSize(3)
        outputVolume(:,:,i) = ~(imfill(logical(volumeShape(:,:,i)),[1 1],4));
    end
    if (all(outputVolume) |  ~any(outputVolume))
        warning('Matlab:surface2volume_nested','Sorry, failed again.');
    end
end

% ovt = outputVolume ;
% outputVolume = logical(imdilate(uint8(outputVolume),strel('ball',1,2))) ;
% outputVolume = double(outputVolume) ;
%@@@ if this assertion never fails, lose the 3 lines above here
% assert(meq(ovt,outputVolume)) ;
%WAITWUT!? strel('ball',1,2) == [1]... so this imdilate (=?=conv) doesnt do
%anything...  easiest 2x speed increase EVER :O
if verboseOutput
    fprintf(outputFID,'done in %g sec\n',toc);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: pointsontriangle                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputPoints = pointsontriangle_nested(triVertices, minElementLength)

DEFAULT_MINELEMENT_LENGTH = 0.25;

% calculate the 3 vectors of the triangle
edges = triVertices([1,1,2],:)-triVertices([2,3,3],:);

% calculate the lengths of these vectors
el = zeros(3,1);
el(1) = norm(edges(1,:));
el(2) = norm(edges(2,:));
el(3) = norm(edges(3,:));

if (~exist('minElementLength','var'))
    minElementLength = DEFAULT_MINELEMENT_LENGTH;
end

% if none of the edges is larger than the minElementLength return the
% inputVertices
if (~(nnz(el > minElementLength)))
    outputPoints = triVertices;
    return
end

longestEdge = max(el);

B = [(0:minElementLength/longestEdge:1) 1];

A = repmat(B,[1 length(B)]) ;
B = reshape(repmat(B,[length(B), 1]),1,[]) ;
C = 1-A-B;

idInTriangle = find(C >= 0);
outputPoints = zeros(length(idInTriangle),3);
outputPoints(:,1) = (triVertices(1,1).*A(idInTriangle) + triVertices(2,1).*B(idInTriangle) + triVertices(3,1).*C(idInTriangle))';
outputPoints(:,2) = (triVertices(1,2).*A(idInTriangle) + triVertices(2,2).*B(idInTriangle) + triVertices(3,2).*C(idInTriangle))';
outputPoints(:,3) = (triVertices(1,3).*A(idInTriangle) + triVertices(2,3).*B(idInTriangle) + triVertices(3,3).*C(idInTriangle))';
end


function [CoMvolume,IAvolume,SVvolume] = InertialAxesVolume_nested(VolumeGrid,PositionGrid)
%INERTIALAXESVOLUME.M determine the principal axes of inertia and center of
% mass of a volume element. Use surface2volume_nested.m to compute the VolumeGrid
% and PositionGrid of a closed surface.
%
% Input:
% VolumeGrid: 0 = no volume voxel 1 = volume voxel [m x n x p matrix]
% PositionGrid: position voxels [1x3 cell -> m x n x p matrix]
%
% Output:
% CoMvolume: center of mass volume object
% IAvolume: principal axes of inertia
% SVvolume: singular values
%
% ORL Nijmegen, september 2012
% Revised (Vectorisation) sept 2016, Max B

% Centre of Mass determination (volume)
gridSize = PositionGrid{1,1}(1,1) - PositionGrid{1,1}(2,1);
I = find(VolumeGrid) ;


CoMxx = PositionGrid{1,1}(I)+ 0.5*gridSize ;
CoMyy = PositionGrid{1,2}(I)+ 0.5*gridSize ;
CoMzz = PositionGrid{1,3}(I)+ 0.5*gridSize ;

CoMvolume = mean([CoMxx CoMyy CoMzz],1);

if nargout > 1
    Ixx = sum((CoMyy-CoMvolume(2)).^2 + (CoMzz-CoMvolume(3)).^2) ;
    Iyy = sum((CoMxx-CoMvolume(1)).^2 + (CoMzz-CoMvolume(3)).^2);
    Izz = sum((CoMxx-CoMvolume(1)).^2 + (CoMyy-CoMvolume(2)).^2);
    
    Ixy = sum((CoMxx-CoMvolume(1)).*(CoMyy-CoMvolume(2)));
    Ixz = sum((CoMxx-CoMvolume(1)).*(CoMzz-CoMvolume(3)));
    Iyz = sum((CoMyy-CoMvolume(2)).*(CoMzz-CoMvolume(3)));
    
    IT = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz]; % moment of inertia tensor
    
    % Principle axes of inertia
    [~,S,IAvolume] = svd(IT);            % axes of inertia
    SVvolume = diag(S);                  % singular values
end
end
