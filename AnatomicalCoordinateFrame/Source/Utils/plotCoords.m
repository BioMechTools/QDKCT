function [ hO ] = plotCoords(varargin)
%PLOTCOORDS Draws a refernce frame in the current axes
%   
%   plotCoords(frame,[origin])
%       frame        -  a 3x3 matrix with each column is a direction
%       origin       -  a vector of length 3 with the coordinates of the origin
%   plotCoords(framestruct)
%       framestruct  -  a struct with 4[3] fields: [origin] and either of [X,Y,Z],
%       [AP PD ML] or [Xdir Ydir Zdir]
%   plotCoords(...,'drawLegend')
%       adds a legend to the current axis with the names of the direction
%   plotCoords(...,'Color',cc)
%       Specify the color of the axes. default is hsv(3) (=rgb)
%   plotCoords(...,'arrowLength')
%       Specify the plotlength of the arrows (20)
%
%   [hO] = plotCoords(...) 
%       hO           -   handles to the shafts of the arrows (which are lines)
%   
%   Example:
%       To plot the unit reference frame (X,Y,Z) in slightly fabulous colors:
%       plotCoords(eye(3),[1 2 3],'Color',cool(3),'drawLegend')
%
%  (V1.01) ORL Nijmegen, Max Bakker 2016

%changelog, added doc, and options for arrowL and legendplotting

axis equal 
axes = eye(3) ;
origin = zeros(1,3) ;
drawLegend = false ;
arrowL = 20 ;

cc = eye(3) ;


if nargin > 0
    ccI = find(strcmp(varargin,'Color')) ;
    if ~isempty(ccI)
        cc = varargin{ccI+1} ;
        varargin([0 1]+ccI)= []; 
    end
    
    axisI = cellfun(@isstruct,varargin) ;
    if nnz(axisI) ==0
        axisI = cellfun(@numel,varargin)==9 ;
    end
    if nnz(axisI)==1
        axes = varargin{axisI} ;
        if(isfield(axes,'origin'))
            origin = axes.origin ;
        elseif isfield(axes,'CM') ;
            origin = axes.CM ;
        end
    end
    varargin(axisI) = [] ;
    
    originI = cellfun(@numel,varargin)==3 ;
    if nnz(originI)==1
        origin = varargin{originI};
        varargin(originI) = []; 
    end
    
    
    legendI = strcmp(varargin,'drawLegend') ;
    if ~isempty(legendI)
        drawLegend = true ; 
        varargin(legendI)= []; 
    end
    
    lengthI = find(strcmp(varargin,'arrowLength')) ;
    if ~isempty(lengthI)
        arrowL = varargin{lengthI+1} ;
        varargin([0 1]+lengthI)= []; 
    end
end

% nDim = numel(origin) ;

if ~isstruct(axes) ;
    axes=  axes' ;
    Naxes.X = axes(1,:) ;
    Naxes.Y = axes(2,:) ;
    Naxes.Z = axes(3,:) ;
    axes = Naxes ;
    clear Naxes
end
fnms = fieldnames(axes) ;
% assert(nDim == 3) ;

%if the struct names are AP, PD and ML, sort them to give the right color
testSystem{1}= {'AP','PD','ML'};
testSystem{2} = {'X','Y','Z'};
testSystem{3} = {'Xdir','Ydir','Zdir'};

for ai = 1 : numel(testSystem) ;
    tisSys = testSystem{ai} ;
    if nnz(any(cell2mat(cellfun(@(s) strcmp(tisSys,s),fnms,'uni',0))))==3
        fnms = tisSys ;
        break ;
    end
end

origin = reshape(origin,1,3) ;

hold on
%number of samples in theta determines arrowhead roundness (now its 10)
nPointsForTheHeads = 10 ;
theta=linspace(0,2*pi,nPointsForTheHeads) ;
%These two might better be parameters maybe
tipLength = arrowL*.1 ;
arrowRadius = tipLength / 2; 
for di = 1 : numel(fnms)
    coord =  reshape(axes.(fnms{di}),1,3) ;
    coords = [origin;origin]+[zeros(1,3);coord]*arrowL ;
    ha(di) = plot3(coords(:,1),coords(:,2),coords(:,3),'Color',cc(di,:),'LineWidth',2) ;
    %arrowheads
    normal = diff(coords,1) ;
    normal = normal / rssq(normal) ;
    v=null(normal) ;
    center = coords(2,:)-tipLength * diff(coords) / rssq(diff(coords)) ;
    points=repmat(center',1,size(theta,2))+arrowRadius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    V = [coords;points'] ;
    F = [ones(1,nPointsForTheHeads+1);2:nPointsForTheHeads+2;[3:nPointsForTheHeads+2,2]]' ;
    patch('Faces',F','Vertices',V,'FaceColor',cc(di,:),'edgeColor',[.4 .4 .4]) ;
end
if drawLegend
    legend(ha,fnms) ;
end
if nargout > 0
    hO = ha ;
end
end

