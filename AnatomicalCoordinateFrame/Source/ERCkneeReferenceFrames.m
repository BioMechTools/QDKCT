function [ refFrames,stlData,DiagInfo ] = ERCkneeReferenceFrames( Femur,Patella,Tibia,axEst)
%ERCKNEEREFERENCEFRAMES Computes the anatomical reference according to the ERC
%standard. A more detailed definition of the standard can be found in the
%help of the respective implementations (see see also below).
%
%  [ refFrames,stlData,DiagInfo ] = ERCKNEEREFERENCEFRAMES( Femur,Patella,Tibia,axEst)
%
%--Input
%   Femur,Patella  -    strings with the paths to the respective stl
%   Tibia               files. OR structs with the fields F and V, holding
%                       the faces and vertices respectively. Any of these
%                       can be empty (== []) ;
%   [axEst]        -    Reference frame estimate. Does not need to be vary 
%                       accurate. If none is supplied the frame will be
%                       estimated depending on the input. Only either 
%                       patella alone or tibia and patella are given the
%                       reference frame is required. axEst is either a
%                       struct with the fields X Y and Z, or a 3x3 matrix
%                       where the first column is X(...). This means that
%                       when the stls are alligned to global space, eye(3)
%                       suffices.
%
%--Output
%   refFrames      -    Struct with the reference frames in the fields f, p
%                       and t. Each field has the fields (X,Y,Z), with the 
%                       coordinates of the respective axis 
%   stlData        -    Surfaces meshes use to determine the frames. 
%   DiagInfo       -    Diagnostics information, use the respective
%                       plotting functions to visualise the steps
%
%  (V1.0) ORL Nijmegen, Max Bakker 2016
%
%  See also Demo ERCrefFrameFemur ERCrefFrameTibia ERCrefFramePatella
%  plotPatellaDiagInfo plotTibiaDiagInfo stlLoad
%  plotFemurDiagInfo eye


%Femur is double (0) when uigetfile is used and cancel is pressed 
%merging file and path results in [0 0], checking for double covers both
if ~exist('Femur','var') || isa(Femur,'double')
    Femur = [];
end
if ~exist('Patella','var')|| isa(Patella,'double')
    Patella = [];
end
if ~exist('Tibia','var')|| isa(Tibia,'double')
    Tibia = [];
end


whichBones = ~cellfun(@isempty,{Femur,Patella,Tibia}) ;
%convert binary length-3 array to numerical value for the switch
whichBones = sum(whichBones .* 2.^(0:2)) ;

%load stuff
if ischar(Femur)
    Femur = stlLoad(Femur) ;
end
if ischar(Patella)
    Patella = stlLoad(Patella) ;
end
if ischar(Tibia)
    Tibia = stlLoad(Tibia) ;
end

%cut the whole femur
if ~isempty(Femur) ;
    [Femur.V,Femur.F,PDest,IAfemur] = wholeFemurCheck_nested(Femur.V,Femur.F) ;
end

%% find the estimate axes

%calculate the estimate frame when it is not supplied
%take into account the different situations
if ~exist('axEst','var')
    
    %2 would be patella only
    %6 would be patella-tibia
    switch whichBones
        
        case 7        %all 3 knee bones
            % add parameterrunners precutting here
            axEst = calculateExternalFrame(Femur.V,Patella.V,Tibia.V) ;
            
        case 5  %femur-tibia
            APest = manualSelectAP_nested(Femur,Tibia,IAfemur) ;
            
        case 4 %tibia only
            [~,IA,~] = computeInertialAxes(Tibia.F,Tibia.V,3) ;
            selection = manualSelectionF((IA'*Tibia.V')','Proximal','Select the Tibial plateau') ;
            PDest = sign(selection*-1+1)*IA(:,1)' ;
            APest = manualSelectAP_nested(Tibia,IA) ;
            
        case 3 %patella-femoral
            APest = manualSelectAP_nested(Femur,Patella,IAfemur) ;
            
        case 1  %Femur
            APest = manualSelectAP_nested(Femur,IAfemur) ;
            
        otherwise
            warning('This combination of bones is not supported') ;
            refFrames= 0 ;
            return
    end
    if whichBones ~=7
        axEst.X = APest ;
        axEst.Y = PDest ;
        axEst.Z = cross(APest,PDest) ;
    end
else
    if ~isstruct(axEst) ;
        aet =axEst ;
        axEst.X = aet(:,1) ;
        axEst.Y = aet(:,2) ;
        axEst.Z = aet(:,3) ;
    end
end

%use a cell to determine nargout in the subfunctions
%   (which is, wether to generate DiagInfo or not)
outArgs = cell(1,1+(nargout>2)) ;
if ~isempty(Femur)
    
    %calculate and apply the reference lenght
    [refL,VcsL,IA] = computeReferenceLength(Femur.V) ;   
    refLinPerc = refL*2 ./ range(VcsL(:,1)) *100  ;
    if sign(dot(axEst.Y,IA(:,1)'))==1
        percRange = [0 refLinPerc] ;        
    else
        percRange = [100-refLinPerc 100] ;
    end
    [Femur.F, VcsL] = cutMeshByPercent(Femur.F,VcsL,percRange,1,1) ;
    %rotate back to original CS
    Femur.V = (IA'\VcsL')' ;
    
    [outArgs{:}] = ERCrefFrameFemur(Femur.F,Femur.V,axEst) ;
    refFrames.f = outArgs{1} ;
    %I could encapsulate this in an if, but this is cheap and inconsequential
    DiagInfo.f = outArgs{end} ;
end
if ~isempty(Patella)
    [outArgs{:}] = ERCrefFramePatella(Patella.F,Patella.V,axEst) ;
    refFrames.p = outArgs{1} ;
    DiagInfo.p = outArgs{end} ;
end

if ~isempty(Tibia)
    [outArgs{:}] = ERCrefFrameTibia(Tibia.F,Tibia.V,axEst) ;
    refFrames.t = outArgs{1} ;
    DiagInfo.t = outArgs{end} ;
end
DiagInfo.version = 1.1; 
if nargout >2
    stlData = [Femur,Patella,Tibia] ;
end
end

function [Vout,Fout,PDest,IA] = wholeFemurCheck_nested(V,F)
% estimate PD from a femurs stl, also cut in half when we have a whole
% femur.
% First determine wether we have a whole femur
% from the Crossectional Area plot. If it has two peaks which are similair
% enough (as measured by 4 loosely tresholded measures). Then look at the
% skewness of the distribution of distances to the center the respective
% end.

%get the inertial axes
% [~,IA,~] = computeInertialAxes(F,V,2) ;
%*get the inertial axes crudely (thnx Marco)
[IA,~] = pca(V) ;

%Find the peaks in crossectiona area plot
[CA,sectionCenters,~] = computeCrossSectionAreaAlongAxis(V,IA,2) ;
[peaks,peaklocs,width,prominence] = findpeaks(CA,sectionCenters,'NPeaks',2,'SortStr','descend') ;

%Qualify the peaks and measure wether we have a full femur
if numel(peaks)==2
    %compare the peaks to determine wether its a whole femur
    peakR = peaks(1)/peaks(2) ;
    peaklocsDiff = abs(peaklocs(1)-peaklocs(2)) ;
    widthR = width(1)/width(2) ;
    promR = prominence(1)/prominence(2) ;
    
    wholeFemur =  promR < 10 && peakR < 3 && peaklocsDiff > 70 && widthR < 4  ;
    disp('Whole Femur detected....cutting') ;
else
    %this likely never happends because CA plot is (neverish) smooth
    wholeFemur = false ;
end


%determine which side is the foot from the distance distribution to the
%means
if wholeFemur
    [~,Vend] = arrayfun(@(i) cutMeshByPercent(F,V,peaklocs(i) + [-1 1] .* width(i) ,1,false,IA'),1:2,'uni',0) ;
    %distance distribution
    endMids = cellfun(@mean,Vend,'uni',0) ;
    
    VdistDistr = cellfun(@(v,mv) rssq(bsxfun(@minus,v,mv),2),Vend,endMids,'uni',0) ;
    skews = cellfun(@skewness,VdistDistr);
    [~,footI] = min(skews) ;
    
    PDest = endMids{(footI==1)+1} - endMids{footI}  ;
    PDest = PDest' / rssq(PDest) ;
    halfPerc = [0 50 100 ];
    %select the half with the foot (use logical as index modifer)
    halfPerc = halfPerc( (1:2) + (peaklocs(footI) > 50) ) ;
    
    
    try
        [Fout,Vout] = cutMeshByPercent(F,V,halfPerc,1,true,IA') ;
    catch
        %because the sruface closing is stupidish, precutting without closing
        %the surface decreases memory usage
        morePerc = (halfPerc) ;
        morePerc(morePerc==50) = 50 - (morePerc(morePerc~=50)-50)/10 ;
        secondHalfPerc = morePerc ;
        if peaklocs(footI) > 50
            secondHalfPerc(halfPerc==50) = 100 - 50/55 * 100 ;
        else
            secondHalfPerc(halfPerc==50) = 50/55 * 100 ;
        end
        [Fout,Vout] = cutMeshByPercent(F,V,morePerc,1,false,IA') ;
        [Fout,Vout] = cutMeshByPercent(Fout,Vout,secondHalfPerc,1,true,IA') ;
    end
    
else
    %PDest is the whole femur's long inertial axis, pointed away from the
    %foot
    PDest = -sign((peaks(1)<50)-.5) .* IA(:,1) ;
    Vout =V;
    Fout = F ;
end

%placeholders :o
% [~,Vwide] = cutMeshByPercent(F,V,peaklocs(footI) + [-1 1] *.25 ,1,false,IA') ;
% APest = condyleContourAPest_nested(Vwide) ;
end

function [APest] = manualSelectAP_nested(varargin)
%gateway to manual selection function
if ~isstruct(varargin{end})
    IA = varargin{end} ;
    varargin  =varargin(1:end-1) ;
else
    [~,IA,~] = computeInertialAxes(varargin{1}.F,varargin{1}.V,3) ;
end
dispV = cell2mat(cellfun(@(s) s.V,varargin,'uni',0)') ;
direction = 3;
selec = manualSelectionF((IA'*dispV')','Posterior','',2,'x',direction) ;
Xest = makehgtform('axisrotate',IA(:,1),selec) * [IA(:,3);0]  ;
%the sign is because we selected posterior, and AP is towards A
APest = -Xest(1:3)' ;
end



% function [AP] = condyleContourAPest_nested(Vwide)
% %Placeholder :O stay tuned
% end
