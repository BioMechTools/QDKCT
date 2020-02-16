function ho = ezplotm(X,varargin)
%Interface to plot3. Basically the same as plot3(X(1,:),X(2,:),X(3,:)), but 
%with automatic detection of array sizes (i.e. allows both ezplotm(X) and
%ezplotm(X') ). 2D is also supported, but 3D is assumed when in doubt.
%   This means that you can not plot 3 2D points with this function.
% Supports X as (with i = 2 or 3) : 
%  -matrix ixN & Nxi
%  -cell arrays per dimension (1xi | ix1 cellarray with 1xN | Nx1 arrays)
%  -AxBxC logical arrays (like a 3d spy)
%
%   see also spy, plot3
%   Max Bakker, November 2016
X = squeeze(X) ;
xs = size(X) ;
if max(xs(xs<4))==2
    nDim = 2;
elseif min(xs)<4
    nDim= 3;
elseif min(xs)>3
    nDim = numel(xs) ;
    [Y(1,:),Y(2,:),Y(3,:)]  = ind2sub(xs,find(X~=0)) ;
    X = Y ;
end
if ~any(cellfun(@(v) isa(v,'char'),varargin))
    varargin(end+1) = {'.'} ;
end


%if a 2D matrix, check which dim holds the coords
if numel(xs)==2
    if xs(2)==nDim
        X =X' ;
    end
end

if numel(X)==nDim && iscell(X)
    reshape(X,[1,nDim]) ;
    if numel(unique(cellfun(@numel,X)))==1 ;
        if size(X{1},1) > size(X{2},2)
            X = cell2mat(X) ;
            X = reshape(X,numel(X)/nDim,nDim)' ;
        else
            X = cell2mat(X) ;
            X = reshape(X,nDim,numel(X)/nDim) ;
        end
    end
end

%one datapoint will be plotted with a dot when no properties are specified
if size(X,2) == 1 && (numel(varargin)==0 || ~any(cellfun(@(v) isa(v,'char'),varargin)))
    varargin(end+1) = {'.'} ;
end

if nDim==2
    h = plot(X(1,:),X(2,:),varargin{:}) ;
else
    h = plot3(X(1,:),X(2,:),X(3,:),varargin{:}) ;
end
if nargout > 0
    ho = h ;
end

end