function plotPlacer(varargin)
%PLOTPLACER, Interface to subplot.
%   Subplot is tedious. plotPlacer only takes arguments the first call, and
%   remembers those for placing all other subplots. Every call without
%   arguments spawns a new subplot, starting anew in the current figure.
%
%   --Input (first call only!)
%      plotPlacer(nPlotx,nPloty).  %plots in this rectangle
%      plotPlacer(nPlot). %plot on a rectangle in 1 figure
%      plotPlacer(...,'vert') draws subsequent axes vertically ('horz' also
%      works. Because of legacy and horizontal and vertical are confusing)

% Max Bakker, September 2014

%Todo?
%add support for non regular placement
%Force overwrite new figures as option?
%add option to supply number of plots and number of figures

persistent nPlotx nPloty fn sp nsp next ;
%place a plot
if nargin == 0
    if(sp == nsp);
        fn = fn + 1 ;
        %call figure here for clf ;
        figure(fn) ;
        clf ;
    end
    %try catches when inplots == 0 and do nothing
    try
        %always call figre for more robust plot placement
        figure(fn);
        sp = next(sp,nPlotx,nPloty) ;
        subplot(nPloty,nPlotx,sp) ;
    catch
    end
else
%process inputs
    if ischar( varargin{end} )
        options = varargin{end} ;
        varargin = varargin(1:end-1) ;
    end
    %initialize
    inSizes = cellfun(@numel,varargin) ;
    
    %non-vararr input
    if all(inSizes==1) && numel(inSizes) < 3
        %add 'horz' to plot them vertically. Now vert also works
        if exist('options','var') && (strcmp(options,'horz') || strcmp(options,'vert') )
            next = @(n,nx,ny) ((n+nx)>nx*ny)*(-nx*ny+1) + n + nx + (n==(nx*ny))*(-nx) ;
            %O yes he just did
        else
            next =  @(n,nx,ny) mod(n, nx * ny)+1 ;
        end
        
        switch numel(varargin)  %~= nargin! bc above the options is removed
            case 1
                nPlotx = ceil(sqrt(varargin{1})) ;
                nPloty = ceil(varargin{1}/nPlotx) ;
            case 2
                nPlotx = varargin{2} ;
                nPloty = varargin{1} ;
        end
    else
        %Hiddem feature. To unify the input with that of parPicker (deprec)
        if exist('options','var') && ~strcmp(options,'horz')
            names= strsplit(options) ;
            varnames = arrayfun(@inputname,1:numel(varargin) ,'uni',0) ;
            
            intmp = cellfun(@(x) strcmp(names,x),varnames,'uni',0) ;
            intmp = cell2mat(reshape(intmp,[numel(intmp),1])) ;
            
            pioritySizes = arrayfun(@(x) inSizes(intmp(:,x)),1:size(intmp,2)) ;
        else
            pioritySizes = inSizes(find(inSizes>1,2,'first'));
        end
        
        next =  @(n,nx,ny) mod(n, nx * ny)+1 ;
        if numel(pioritySizes)>0, nPlotx = pioritySizes(1) ; else  nPlotx = 1 ; end
        if numel(pioritySizes)>1, nPloty = pioritySizes(2) ; else  nPloty = 1 ; end
    end
    fig = gcf ;
    fn = get(fig,'Number') ;
    if ischar(fn)
        %Old matlab compatibility
        fn = fig ;
    end
    fn = fn - 1 ;%-1 because we want a call to figure first time
    nsp = nPlotx * nPloty;
    sp = nsp ;
end
end