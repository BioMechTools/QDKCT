function [refFrame,DiagInfo] = ERCrefFrameTibia(F, V,axEst,Param)
%ERCREFFRAMETIBIA(F,V,axEst,Param) computes the ERC standard anatomical 
%   reference frame for a proximal tibia. 
%
%   The ERC standard reference frame is defined according to Miranda e.a.'s
%   definition. The tibial plateau's inertial axes are sorted based on
%   their moment of inertia. From heigh to low moment of inertia they give: 
%       PD(Y): Largest plateau moment
%       AP(X): Left as an exercise to the reader
%       ML(Z): Smallest plateau moment
%--Input
%   F        -    Faces of the tibia
%   V        -    Vertices of the tibia
%   axEst    -    Axis estimate
%                    Used to determining the signs of the axes
%   Param    -    Parameter structure with (any of the) fields( = default):
%                      Param.sectionStepSize = 2 ;
%                           Step size in percentage for determining the
%                           crossection area
%                      Param.IAgridSize = [0.5 1] ;
%                          point-mass-grid-size for determining the
%                          inertial axes
%
%--Output
%   refFrame  -   Struct with three field (X,Y,Z), having the coordinates 
%                 of the respective axis 
%   DiagInfo -    Diagnostics information, use plotTibiaDiagInfo to show
%                     the different steps.
%
%
%  ORL Nijmegen, WJ Zevenbergen September 2012
%  (V1.0) ORL Nijmegen, Max Bakker 2016
%
% Reference:
%   -Miranda, D. L., M. J. Rainbow, et al. (2010). "Automatic determination of
%    anatomical coordinate systems for three-dimensional bone models of the
%    isolated human knee." J Biomech 43(8): 1623-1626.

%% Initialising
try
    P_crossectStepSize = 2 ;
    P_IAgridSize = [0.5 1] ;
    if nargin > 3
        parfields  = fieldnames(Param) ;
        for pf = parfields'
            switch pf{1}
                case 'sectionStepSize'
                    P_crossectStepSize = Param.sectionStepSize ;
                    
                case 'IAgridSize'
                    P_IAgridSize = Param.IAgridSize ;
            end
        end
    end
    if numel(P_IAgridSize)==1
        P_IAgridSize = repmat(P_IAgridSize,[2,1]) ;
    end
    
    
    if nargout > 1
        DiagInfo.param.sectionStepSize = P_crossectStepSize ;
        DiagInfo.param.IAgridSize = P_IAgridSize ;
    end
    
    %% the algorithm
    % determine center of mass and inertial axes
    [CoM, IA,SV] = computeInertialAxes(F,V,P_IAgridSize(1));
    % calculate crossectional area (CA) along the long axis
    [CA,sectionCenters,CMonCAaxis] = computeCrossSectionAreaAlongAxis(V,IA,P_crossectStepSize) ;
    
    % We seperate the tibial plateau from the tibia at the 'widest' (largest CA) point
    [~,iCA1] = max(CA);         % maximal cross-sectional area
    CA1_perc = sectionCenters(iCA1);       % percentage max cross-sectional area
    
    % Selecting tibia plateau
    plateauRange_perc = [0 CA1_perc 100] ;
    plateauRange_perc  = plateauRange_perc([1 2] + (CMonCAaxis < CA1_perc) ) ;
    %compute inertial axes of the plateau
    [CoM_plateau, IA_plateau, SV_plateau,V_plateau,F_plateau] = computeInertialAxes(F,V,P_IAgridSize(2),plateauRange_perc,IA) ;
    
    %the external frame determines the sign of the ref frame
    %the direction is determined by the plateau's IA
    %note that IA directions are sorted based on the moment
    Miranda_origin = CoM_plateau ;
    Miranda_APaxis = sign(dot(axEst.X,IA_plateau(:,2))) * IA_plateau(:,2)';
    Miranda_PDaxis = sign(dot(axEst.Y,IA_plateau(:,3))) * IA_plateau(:,3)';
    Miranda_MLaxis = sign(dot(axEst.Z,IA_plateau(:,1))) * IA_plateau(:,1)';
    
    %% Output
    % output, and chose the positive direction according to the input axes
    refFrame.X = Miranda_APaxis;
    refFrame.Y = Miranda_PDaxis;
    refFrame.Z = Miranda_MLaxis;
    refFrame.origin = Miranda_origin;
    
    % populate the diagnostics information structure
    if nargout > 1
        DiagInfo.IA.Tibia.IA = IA ;
        DiagInfo.IA.Tibia.SV = SV ;
        DiagInfo.IA.Tibia.CM = CoM ;
        
        DiagInfo.crossSect.percCM = CMonCAaxis ;
        DiagInfo.crossSect.section = sectionCenters ;
        DiagInfo.crossSect.data = CA ;
        DiagInfo.crossSect.i1 = iCA1 ;
        %     DiagInfo.crossSect.Usage = 'Use plot((...).section(i1),(...).Data(i1)) to add the dot of P1';
        
        DiagInfo.sections.Plateau.F = F_plateau ;
        DiagInfo.sections.Plateau.V = V_plateau ;        % vertices tibia plateau
        
        DiagInfo.IA.Plateau.IA = IA_plateau ;
        DiagInfo.IA.Plateau.SV = SV_plateau ;
        DiagInfo.IA.Plateau.CM = CoM_plateau ;
        
        Vrot = (IA' * V')' ;
        DiagInfo.sections.Plateau.cutPlane  = reshape((IA' \ [CA1_perc/100*range(Vrot(:,1))+min(Vrot(:,1)) 0 0; 0 1 0; 0 0 1  ] ),[1,9]) ;
        
        DiagInfo.refFrame = refFrame ;
    end
    DiagInfo.failed = 0;
catch emAll
    % For instance to short diaphysis causes the method to crash
    %this block makes sure the output still is somewhat sensible
    emAll.getReport
    
    DiagInfo.failed = 1;
    DiagInfo.errorMessage = emAll ;
    
    %phony data to reveal it goes wrong w/o crashing
    dummyCoords.X = [1 0 0] ;
    dummyCoords.Y = [-1 0  0] ;
    dummyCoords.Z = [0 1 1] ;
    dummyCoords.origin = [ 0 0 0 ] ;
    DiagInfo.refFrame = dummyCoords ;
end
end
