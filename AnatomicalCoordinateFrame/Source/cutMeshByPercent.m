function [F_sel, V_sel] = cutMeshByPercent(F,V,section_perc,dim,closeSurface,rotM)
%CUTMESHBYPERCENT select a part of an stl-file.
%  [F_sel, V_sel] = cutMeshByPercent(F,V,section_perc,dim,niceEdges,rotM)
%   Cuts a specific part from a mesh, within a certain range (in
%   percentage) and a certain diminsion. Using state-of-the-art-technology
%   vertices can be moved to result in a nice cutting-edge. A rotation matrix
%   may be supplied (with directions in columns) to specify a non-primary
%   direction.
%
%---Input
%     F              -    Faces (N by 3 integer
%     V              -    Vertices (N by 3)
%     section_perc   -    Section range in percentage (e.g. [80 100]) 
%     dim            -    Dimension in which to cut ([1,2,3])
%     [closeSurface] -  Flag wether to close the surface(logical [false]) 
%          Using the rule: out-of-range-vertices in faces which also have
%          in-range-vertices are moved towards the mean connected
%          in-range-vertex until they are on the edge.
%     [rotM]       -    Rotation matrix defining the coordinate system
%          When rotM is supplied, dim is about rotM's columns
%
% Output:
%     F_sel          -     Faces of the section
%     V_sel          -     Vertices of the section
%
%  ORL Nijmegen, WJ Zevenbergen May 2012
%  (V1.0) ORL Nijmegen, Max Bakker 2016
%     
%   see also : computeInertialAxes
if min(section_perc) >=100 || max(section_perc)<=0;
    error('STL section input percentage is impossible') ;
end

if nargin < 5 || isempty(closeSurface)
    closeSurface = false ;
end

% unique vertices of the object
[V, ~, indexn] = unique(V, 'rows'); F = indexn(F); clear indexn
if nargin > 5
    V = (rotM * V')';
end

% vertice range
minV =  min(V(:,dim));

secInCoords =  section_perc/100*( max(V(:,dim))-min(V(:,dim))) ;

% index selected vertices
nV_sel = find(V(:,dim) <= minV + secInCoords(2) & V(:,dim) >= minV + secInCoords(1));


% Faces
nV_sel_F_sel_check1 = zeros(size(F));      % empty matrix
for k = 1:length(nV_sel)
    nV_sel_F_sel = (F == nV_sel(k));       % find selected vertice #k in face matrix
    nV_sel_F_sel_check1 = nV_sel_F_sel_check1 + nV_sel_F_sel; % selected vertices in face matrix
end
targetHeights(1) = max(V(nV_sel,dim) ) ;
targetHeights(2) = min(V(nV_sel,dim) ) ;
NvertPerFace = sum(nV_sel_F_sel_check1,2) ;

nFace_sel = NvertPerFace> 2; % true face is selected (is different with nice edges)

%move vertices to the plane when they are on the wrong side AND are in
%faces which have vertices on the right side
if closeSurface
    correctTheseFaces = NvertPerFace>0 & NvertPerFace <3 ;
    if nnz(correctTheseFaces)>0
        FtoCorrect = F(correctTheseFaces,:) ;
        correctTheseVertices = unique(FtoCorrect .* ~nV_sel_F_sel_check1(correctTheseFaces,:),'sorted') ;
        
        if correctTheseVertices(1) == 0
            correctTheseVertices = correctTheseVertices(2:end) ;
        end
        for vi =correctTheseVertices'
            refVerts = unique(FtoCorrect(any(FtoCorrect==vi,2),:)) ;
            refVerts = setdiff(refVerts,correctTheseVertices) ;
            if numel(refVerts)>1
                useTheseV = ~bsxfun(@eq,bsxfun(@le, V(refVerts,dim), targetHeights ),V(vi,dim) < targetHeights ) ;
                assert(nnz(useTheseV)>0) ;
            else
                useTheseV = [1 1];
            end
            
            goodVertRef = mean(V(refVerts(useTheseV(:,1)),:),1) ;
            badVert = V(vi,:) ;
            %move bad vertices towards the good vertex (or mean good vertex) in
            %order to end them up at the boundary
            a = (targetHeights(1)- badVert(:,dim)) ./ (goodVertRef(:,dim) - badVert(:,dim)) ;
            if any(abs(a)>1) || isnan(a)
                goodVertRef = mean(V(refVerts(useTheseV(:,2)),:),1) ;
                a= (targetHeights(2)- badVert(:,dim)) ./ (goodVertRef(:,dim) - badVert(:,dim)) ;
            end
            V(vi,:) = bsxfun(@times,(1-a), badVert) + bsxfun(@times,a,goodVertRef) ;
            %                     assert(~any(any(isnan(V))))
        end
        nV_sel = unique([nV_sel ;correctTheseVertices]) ;
        nFace_sel = NvertPerFace> 0; % true face is selected
    end
end


V_sel = V(nV_sel,:) ;


F_sel = [];
n = 1;
for i = 1:size(F,1)
    if nFace_sel(i) == 1;
        c1 = find(nV_sel == F(i,1)); % 1st vertice face
        c2 = find(nV_sel == F(i,2)); % 2nd vertice face
        c3 = find(nV_sel == F(i,3)); % 3rd vertice face
        F_sel(n,:) = [c1, c2, c3];   % selected face
        n = n + 1;
    end
end


% uniqueify the vertices. ununiqueness might be caused by the moving of vertices
[V_sel, ~, indexn] = unique(V_sel, 'rows'); F_sel = indexn(F_sel);
%when two vertices get merged, a face becomes degenerate (one vertex twice)
F_sel = F_sel(~(F_sel(:,1)==F_sel(:,2)|F_sel(:,2)==F_sel(:,3)|F_sel(:,3)==F_sel(:,1)),:) ;


if(closeSurface)
    [F_sel, V_sel] = CloseSurface_nested(F_sel,V_sel,dim,.5) ;
end

%potentially rotate back to the original reference frame
if nargin > 5
    V_sel = (rotM \ V_sel')' ;
end
end


function [F_enclosed, V_enclosed] = CloseSurface_nested(F,V,dim,closeRange)
%CLOSESURFACE.M function to close a surface which is previously selected
% using cutMeshByPercent.m. Closing this surface is necessary to convert
% the selected surface to a volume object using (surface2volume_nested.m)
%
% Input:
% F - Faces object: contains the vertex lists defining each triangle face [n x 3]
% V - Vertices object: contains the vertices for all triangles [3*n x 3]
% ax_section - Section axis (1 = x, 2 = y or 3 = z)
%
% Output:
% F_enclosed - Faces closed surface
% V_enclosed - Vertices closed surface unique
%
% ORL Nijmegen, WJ Zevenbergen, September 2012
% made a really closed surface on sept '16

% determine plane to close the surface
ax_plane = setdiff(1:3,dim) ;

% find edge vertices
%added by max;
%find edge poins and edges
allCombs =[1 2; 1 2 ; 3 3]' ;
%construct the list of edges (size will be *3 number of faces by 2)
FindMatV = reshape(F(:,allCombs),[numel(F),2]) ;
%count the edges using sparse to construct the edgecountermatrix
edgeOccurCounterV = sparse(FindMatV(:,1),FindMatV(:,2),ones(size(FindMatV,1),1),max(F(:)),max(F(:))) ;
edgeOccurCounterV = edgeOccurCounterV + edgeOccurCounterV' ;
BorderVertices = any(edgeOccurCounterV==1,2) ;

% borderEdges = cell(1,2); [borderEdges{:}] = find(edgeOccurCounterV == 1) ;
% borderEdges = unique(sort([borderEdges{:}],2),'rows') ;

% minimal and maximal vertex values
minVPos = min(V,[],1);
maxVPos = max(V,[],1);

% Select minimal vertex values use to close the surface with a DelaunayTriangulation
range_minV = [minVPos(dim),minVPos(dim) + closeRange];
nV_selmin = V(:,dim) >= range_minV(1) & V(:,dim) <= range_minV(2);
selTheseMin =nV_selmin & BorderVertices ; 
V_selmin = V(selTheseMin,:);
% DTmin = DelaunayTri(V_selmin(:,ax_plane(1)),V_selmin(:,ax_plane(2)));
% FacesMin = DTmin.Triangulation;
triangulationMin = delaunayTriangulation(V_selmin(:,ax_plane));%,borderEdges) ;
FacesMin = triangulationMin.ConnectivityList ;

F_enclosed_min = [F;(FacesMin + size(V,1))]; % Add faces to input faces
V_enclosed_min = [V;V_selmin];               % Add vertices to input vertices

% Select maximal vertex values use to close the surface with a DelaunayTriangulation
range_maxV = [maxVPos(dim)-closeRange,maxVPos(dim)];
nV_selmax = V(:,dim) >= range_maxV(1) & V(:,dim) <= range_maxV(2);
V_selmax = V(nV_selmax & BorderVertices,:);


triangulationMax = delaunayTriangulation(V_selmax(:,ax_plane)); %,borderEdges) ;
FacesMax = triangulationMax.ConnectivityList ;

F_enclosed_minmax = [F_enclosed_min; (FacesMax + size(V_enclosed_min,1))];
V_enclosed_minmax = [V_enclosed_min; V_selmax];

% remove duplicate vertices
[V_enclosed, ~, indexn] = unique(V_enclosed_minmax, 'rows');
F_enclosed = indexn(F_enclosed_minmax);

F_enclosed = TrimEdges_nested(F_enclosed) ;
end




function [ F ] = TrimEdges_nested( F )
%TRIMEDGES Summary of this function goes here
%   Detailed explanation goes here
%TODO:


%remove faces outside of the original shape by iteratively removing
%borderfaces (faces which have an edge which occurs in only 1 face)
%note that this will stop at the original edge, bc the original shape
%is intact and closed.
%also note that this will remove all faces of a mesh with a hole

allCombs =[1 2; 1 2 ; 3 3]' ;

qi = 1 ;
while(true)
    %THIS IS quite stupid (in the sense that it could be smarter, but maybe this is fast)
    %currently, check each face every time
    %we should (maybe?!) check faces neighbouring a recentely deleted face
    %in all but the first iteration.
    findMat = reshape(F(:,allCombs),[numel(F),2]) ;
    %count the edges using sparse to construct the edgecountermatrix
    edgeOccurCounter = sparse(findMat(:,1),findMat(:,2),ones(size(findMat,1),1),max(F(:)),max(F(:))) ;
    edgeOccurCounter = edgeOccurCounter + edgeOccurCounter' ;
    [ii2,jj2] = find(edgeOccurCounter == 1) ;
    borderEdges = [ii2,jj2] ;
    %remove duplicates (directional insensitive)
    borderEdges = unique(sort(borderEdges,2),'rows');
    fwe1 = any(bsxfun(@eq,F,reshape(borderEdges(:,1),[1 1 size(borderEdges,1)])),2) ;
    fwe2 = any(bsxfun(@eq,F,reshape(borderEdges(:,2),[1 1 size(borderEdges,1)])),2) ;
    %Tthese should be removed
    stupidFaces = any(fwe1 & fwe2,3) ;     
    
    if nnz(stupidFaces)==0 %|| qi > 4
        break
    end
    
    F = F(~stupidFaces,:) ;
    qi = qi +1 ;
end
if isempty(F)
    error('Every Face is removed during trimming of the edges. Is there a hole somewhere?') ;
end


end