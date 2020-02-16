function [femCoords, DiagInfo] = ERCrefFrameFemur(F, V, axEst, Param)
%ERCREFFRAMEFEMUR computes the ERC standard anotomical reference frame 
%   for a distal femur.
%
%   The ERC standard reference frame is defined according to Lianne's
%   definition:
%       ML(Z): Line connecting the cylinders fitting the articulating
%              surfaces
%       PD(Y): Diaphysial inertial axis.*
%       AP(X): Perpendicular to both
%   *Almost, the closest vertex perpendicular to ML is used.
%   
%   [femCoords, DiagInfo] = ERCREFFRAMEFEMUR(F,V,axEst,Param)
%
%--Input
%   F        -    Faces of the femur
%   V        -    Vertices of the femur
%   axEst    -    Axis estimate
%                   Used for determining the signs of the axes and choosing
%                   the directions in which to find intersections
%   Param    -    Parameter structure with (any of the) fields( = default):
%                     Param.infPoleFindScaleF = 1.5 ;
%                          Scaling factor for the inferior pole finding
%                     Param.anteriorPercentage = 30 ;
%                          What anterior percentage is used to find AP
%                     Param.IAgridSize = 0.5 ;
%                          point-mass-grid-size for determining the
%                          inertial axes
%
%--Output
%   refFrame  -   Struct with three field (X,Y,Z), having the coordinates of the
%                 respective axis 
%   DiagInfo -    Diagnostics information, use plotFemurDiagInfo to show
%                 the different steps.
%
%
% ORL Nijmegen, WJ Zevenbergen September 2012
%  (V1.0) ORL Nijmegen, Max Bakker 2016
%
% References:
%   -Miranda, D. L., M. J. Rainbow, et al. (2010). "Automatic determination
%    of anatomical coordinate systems for three-dimensional bone models of
%    the isolated human knee." J Biomech 43(8): 1623-1626.
%   -Yin L, Chen K, Guo L, Cheng L, Wang F, et al. (2015) Identifying the 
%    Functional Flexion-extension Axis of the Knee: An In-Vivo Kinematics
%    Study. PLOS ONE 10(6): e0128877. doi: 10.1371/journal.pone.0128877
%
% see also ERCrefFrameTibia ERCrefFramePatella plotFemurDiagInfo

%% Initialising
try
    %set default parameters, overwrite when prsent in the param struct
    P_sectionStepSize = 2;
    P_FirstFitUseEvery = 1 ;
    
    %with in comments tested values which are worse
    P_vertNormTol = 0.6 ;%[ 0.55 0.6 0.65] ;
    P_backDistFacTol = .35 ;%[0.3 0.35 0.4] ;%1/5 ;
    P_midDistFacTol = .24 ;%[.2 .24 .28] %1/6 ;
    P_axDistPostFacTol = .725 ;%%[0.7 0.725 0.75] ;%1/1.7;
    
    P_tolp = 0.1;
    P_tolg = 0.1;
    P_IAgridSize = 0.5 ;
    
    %constants
    C_p3Const = 1.3 ;
    C_p2Const = 0.5 ;
    
    if nargin > 3 && ~isempty(Param)
        parfields  = fieldnames(Param) ;
        for pf = parfields'
            switch pf{1}
                case 'sectionStepSize'
                    P_sectionStepSize = Param.sectionStepSize ;
%                 case 'p3Const'
%                     C_p3Const = Param.p3Const ;
%                 case 'p2Const'
%                     C_p2Const = Param.p2Const ;
                case 'FirstFitUseEvery'
                    P_FirstFitUseEvery = Param.FirstFitUseEvery ;
                case 'vertNormTol'
                    P_vertNormTol = Param.vertNormTol ;
                case 'backDistFacTol'
                    P_backDistFacTol = Param.backDistFacTol ;
                case 'midDistFacTol'
                    P_midDistFacTol = Param.midDistFacTol ;
                case 'axDistPostFacTol'
                    P_axDistPostFacTol = Param.axDistPostFacTol ;
                case 'tolp'
                    P_tolp = Param.tolp ;
                case 'tolg'
                    P_tolg = Param.tolg ;
                case 'IAgridSize'
                    P_IAgridSize = Param.IAgridSize ;
                otherwise
                    warning(['Unkown parameter ' pf{1}]) ;
                    
            end
        end
    end
    %enable different tollerences for the three (1+1+2) different fits
    if numel(P_tolp)==1
        P_tolp = repmat(P_tolp,3,1) ;
    end
    if numel(P_tolg)==1
        P_tolg = repmat(P_tolg,3,1) ;
    end
    
    if numel(P_IAgridSize)==1
        P_IAgridSize = repmat(P_IAgridSize,2,1) ;
    end
    if nargout > 1
        DiagInfo.param.sectionStepSize = P_sectionStepSize ;
        DiagInfo.param.p3Const = C_p3Const ;
        DiagInfo.param.p2Const = C_p2Const ;
        DiagInfo.param.FirstFitUseEvery = P_FirstFitUseEvery ;
        DiagInfo.param.vertNormTol = P_vertNormTol ;
        DiagInfo.param.backDistFacTol = P_backDistFacTol ;
        DiagInfo.param.midDistFacTol = P_midDistFacTol ;
        DiagInfo.param.axDistPostFacTol = P_axDistPostFacTol ;
        DiagInfo.param.tolp = P_tolp ;
        DiagInfo.param.tolg = P_tolg ;
        DiagInfo.param.IAgridSize = P_IAgridSize ;
    end
    
    %% Intertial Axes and Centre of Mass femura
    [CM,IA,SV] = computeInertialAxes(F,V,P_IAgridSize(1)) ;
    
    if nargout > 1 
        DiagInfo.IA.Femur.IA = IA ;
        DiagInfo.IA.Femur.CM = CM ;
        DiagInfo.IA.Femur.SV = SV ;
        
        DiagInfo.param.sectionStepSize = P_sectionStepSize ;
        DiagInfo.param.FirstFitUseEvery = P_FirstFitUseEvery ;
        DiagInfo.param.vertNormTol = P_vertNormTol ;
        DiagInfo.param.backDistFacTol = P_backDistFacTol ;
        DiagInfo.param.midDistFacTol = P_midDistFacTol ;
        DiagInfo.param.axDistPostFacTol = P_axDistPostFacTol ;
        DiagInfo.param.tolp = P_tolp ;
        DiagInfo.param.tolg = P_tolg ;
        DiagInfo.param.IAgridSize = P_IAgridSize ;
        
    
        
    end
    %% Cross sectional area along the smallest inertial axis of the femur (axial length)
    % Normal to the axial length of the femur
    % Rotation of the vertices with respect to the orientation of the inertial
    % axes of the femura
    R_CS1 = eye(3)*IA';            % Rotation matrix inertial axes femura
    V_CS1 = (R_CS1*V')';    % Re-orientation of the femoral vertices
    
    % Calculating axial cross sectional area
    [CA,sectionCenters,CMonCAaxis] = computeCrossSectionAreaAlongAxis(V,IA,P_sectionStepSize) ;
    
    % Locations of the cross-sectional area for the isolation of the condyles
    [CA1, iCA1] = max(CA);
    CA1_perc = sectionCenters(iCA1);
    
    %if statement to handle the different signs (positive directions) of the IA axes
    if CA1_perc <= 50;
        iCA2 = find(CA >= C_p2Const*CA1,1,'last'); % 1/2 max cross-sectional area
        CA2_perc = sectionCenters(iCA2);       % percentage 1/2 max cross-sectional area
        CA3_perc = C_p3Const*CA2_perc;           % percentage 1/2 max cross-sectional area + 30% of %1/2 max CSA
        iCA3 = iCA2 + round((C_p3Const-1)/P_sectionStepSize*CA2_perc);
    else
        iCA2 = find(CA >= C_p2Const*CA1,1,'first');
        CA2_perc = sectionCenters(iCA2);
        CA3_perc = 100-C_p3Const*(100-CA2_perc);
        iCA3 = iCA2 - round((C_p3Const-1)/P_sectionStepSize*(100-CA2_perc));
    end
    
    if nargout > 1
        DiagInfo.crossSect.percCM = CMonCAaxis ;
        DiagInfo.crossSect.section = sectionCenters ;
        DiagInfo.crossSect.data = CA ;
        DiagInfo.crossSect.i1 = iCA1 ;
        DiagInfo.crossSect.i2 = iCA2 ;
        DiagInfo.crossSect.i3 = iCA3 ;
        %         diagInfo.crossSect.Usage = 'Use plot((...).section(i1),(...).data(i1)) to add the dot of P1';
    end
    %% pt1
    
    % diaphysis interial axes determination
    diaphysisRange_perc = [0 CA3_perc 100] ;
    %this is a conditional statement that removes either the 0 or the 100
    %   making sure not to include the condyles (around Point 1)
    diaphysisRange_perc = diaphysisRange_perc([1,2]+(CA1_perc <= 50)) ;
    [CoM_Diaphysis,IA_Diaphysis,SV_Diaphysis,V_Diaphysis,F_Diaphysis] = computeInertialAxes(F,V,P_IAgridSize(2) ,diaphysisRange_perc,IA) ;
    
    % Determine position pt1, using axEst to find the sign towards distal
    pt1 = computeIntersectionLineSurface(V, sign(dot(-axEst.Y,IA_Diaphysis(:,1))) * IA_Diaphysis(:,1),CoM_Diaphysis,0) ;
    if nargout > 1
        DiagInfo.IA.diaphysis.CM = CoM_Diaphysis;
        DiagInfo.IA.diaphysis.IA = IA_Diaphysis ;
        DiagInfo.IA.diaphysis.SV = SV_Diaphysis' ;
        DiagInfo.sections.diaphysis.F = F_Diaphysis ;
        DiagInfo.sections.diaphysis.V = V_Diaphysis ;
        DiagInfo.points12.pt1 = pt1 ;
        
    end
    
    %% pt2
    % Axial plane famur at 2nd location cross-sectional area (1/2 max CSA)
    [~,V_edge_pt2] = cutMeshByPercent(F,V_CS1,[CA2_perc-0.5*P_sectionStepSize CA2_perc+0.5*P_sectionStepSize],1);
    % Center of a box bounding (parallel to diaphyal IA axes) the femur at 2nd location
    boxMiddle_pt2_CS1 = min(V_edge_pt2) + range(V_edge_pt2)/2 ;
    boxMiddle_pt2 = (R_CS1\boxMiddle_pt2_CS1')';
    
    % Pt2 determination
    % Intersection of the femur's inertial axis pointing in AP direction, originated
    % at a box bounding the femur, through the posterior femur
    pt2 = computeIntersectionLineSurface(V, sign(dot(IA(:,3),-axEst.X'))*IA(:,3),boxMiddle_pt2,0);
    if nargout > 1
        DiagInfo.points12.pt2 = pt2 ;
        DiagInfo.points12.pt2_planecenter = boxMiddle_pt2 ;
        
        DiagInfo.points12.pt2_edge= (R_CS1 \ V_edge_pt2')' ;
        DiagInfo.points12.pt2_plane = reshape((R_CS1\[CA2_perc/100*range(V_CS1(:,1))+min(V_CS1(:,1)) 0 0; 0 1 0; 0 0 1  ]'),[1,9]) ;
    end
    %% Isolate condyles (first iteration)
    % Plane's first iteration
    pt1pt2axis = (pt2 - pt1) ;
    pt1pt2axis  = pt1pt2axis / rssq(pt1pt2axis ) ; % Normalize
    MLaxisEstFromIA = IA(:,2)';              % ML axis, inertial axis femur (2)
    % normal vector perpendicular to pt1pt2 axis and ML axis
    APest_iter1 = cross(pt1pt2axis,MLaxisEstFromIA) ;
    APest_iter1 = APest_iter1  / rssq(APest_iter1 ) ;
    %allign with external axis to ensure the condyles are in the first half
    %(aka the 0 in cutMeshByPercent's input is correct)
    APest_iter1 = sign(dot(APest_iter1,axEst.X)) .* APest_iter1 ;
    
    % Rotation to new CS (condyles)
    R_CS2 = eye(3)*[MLaxisEstFromIA' pt1pt2axis' APest_iter1']'; % Rotation matrix from femur CS to condyle CS
    V_CS2 = (R_CS2*V')';                     % Vertices in Condyle CS
    pt1_CS2 = (R_CS2*pt1')';
    
    % Isolating condyles
    pt1pt2_perc_CS2 = (pt1_CS2(3)-min(V_CS2(:,3)))*100/range(V_CS2(:,3));
    [F_condyles_it1, V_condyles_CS2] = cutMeshByPercent(F,V_CS2,[0 pt1pt2_perc_CS2],3);
    % Project back to original coordinate system
    V_condyles_it1 = (R_CS2\V_condyles_CS2')';
    
    if nargout > 1
        DiagInfo.sections.condyles1.cutPlane=reshape((R_CS2\[0 0 pt1pt2_perc_CS2/100*range(V_CS2(:,3))+min(V_CS2(:,3)); 1 0 0; 0 1 0 ]'),[1,9]) ;
        DiagInfo.sections.condyles1.F = F_condyles_it1 ;
        DiagInfo.sections.condyles1.V = V_condyles_it1 ;
        
        DiagInfo.points12.pt2_centerplane = [pt1pt2axis;MLaxisEstFromIA;APest_iter1]'  ;
        DiagInfo.rotationMatrices.CS2 = R_CS2 ;
    end
    %% Cylinder fitting first iteration on isolated condyles
    
    % Determine initial conditions of the cylinder fit
    CM_condyles = mean(V_condyles_it1) ;
    radiusGuess = mean(range(V_condyles_CS2(:,[2,3])))/2 ;
    
    [pointCondyles_i1, axisCondyles_i1, radiusCondyles_i1, ...
        CondyleFit_i1.d, CondyleFit_i1.sigmah, CondyleFit_i1.conv, ...
        CondyleFit_i1.Vx0n, CondyleFit_i1.Van, CondyleFit_i1.urn, ...
        CondyleFit_i1.GNlog, CondyleFit_i1.a, CondyleFit_i1.R0, ...
        CondyleFit_i1.R] = lscylinder(V_condyles_it1(1:P_FirstFitUseEvery:end,:), CM_condyles', MLaxisEstFromIA, radiusGuess, P_tolp(1), P_tolg(1)) ;
    
    if nargout > 1
        DiagInfo.fitCylinders.Miranda1.center = pointCondyles_i1 ;
        DiagInfo.fitCylinders.Miranda1.axis = axisCondyles_i1 ;
        DiagInfo.fitCylinders.Miranda1.radius = radiusCondyles_i1 ;
        DiagInfo.fitCylinders.Miranda1.fitInfo = CondyleFit_i1 ;
    end
    %% Isolate condyles 2nd iteration, improve condyle cutting with first cyl's axis
    % And Cylinder fitting second iteration isolation condyles (Miranda)
    
    % Plane's 2nd iteration
    APest_iter2 = cross(pt1pt2axis,axisCondyles_i1) ;
    APest_iter2 = APest_iter2  / rssq(APest_iter2 ) ;
    %allign with external axis to ensure the condyles are in the first half
    %(aka the 0 in cutMeshByPercent's input is correct)
    APest_iter2 = sign(dot(APest_iter2,axEst.X)) .* APest_iter2 ;
    
    % Rotation to new CS (condyles, but better)
    R_CS3 = eye(3)*[axisCondyles_i1 pt1pt2axis' APest_iter2']'; % Rotation matrix from femur CS to condyle CS
    V_CS3 = (R_CS3*V')';                           % Vertices in Condyle CS
    pt1_CS3 = (R_CS3*pt1')';
    
    % The new plane crosses PT1 (and 2 too) in CS3 and is aligned with axes
    pt1pt2_perc_CS3 = (pt1_CS3(3)-min(V_CS3(:,3)))*100/range(V_CS3(:,3));
    [F_Condyles_it2, V_Condyles2_CS3] = cutMeshByPercent(F,V_CS3,[0 pt1pt2_perc_CS3] ,3);
    
    % Project back to original coordinate system
    V_Condyles_it2 = (R_CS3 \V_Condyles2_CS3')';   % selected condyles in femur CS
    
    %fit a cylinder
    [pointCondyles_i2, axisCondyles_i2, radiusCondyles_i2, ...
        CondyleFit_i2.d, CondyleFit_i2.sigmah, CondyleFit_i2.conv, ...
        CondyleFit_i2.Vx0n, CondyleFit_i2.Van, CondyleFit_i2.urn, ...
        CondyleFit_i2.GNlog, CondyleFit_i2.a, CondyleFit_i2.R0, ...
        CondyleFit_i2.R] = lscylinder(V_Condyles_it2, pointCondyles_i1, axisCondyles_i1, radiusCondyles_i1, P_tolp(2), P_tolg(2));
    strut_Cylinder.pointCondyles_i2=pointCondyles_i2;
    strut_Cylinder.axisCondyles_i2=axisCondyles_i2;
    strut_Cylinder.radiusCondyles_i2=radiusCondyles_i2;
%     save('D:\project\4DTTTG\data\001\strut_Cylinder', 'strut_Cylinder'); 
    if nargout > 1
        %export percplane_condyles2 as a plane according to geom3d's definition
        DiagInfo.sections.condyles2.cutPlane=reshape((R_CS3\[0 0 pt1pt2_perc_CS3/100*range(V_CS3(:,3))+min(V_CS3(:,3)); 1 0 0; 0 1 0 ]'),[1,9]) ;
        DiagInfo.sections.condyles2.V = V_Condyles_it2 ;
        DiagInfo.sections.condyles2.F = F_Condyles_it2 ;
        
        DiagInfo.fitCylinders.Miranda2.center = pointCondyles_i2 ;
        DiagInfo.fitCylinders.Miranda2.axis = axisCondyles_i2 ;
        DiagInfo.fitCylinders.Miranda2.radius = radiusCondyles_i2 ;
        DiagInfo.fitCylinders.Miranda2.fitInfo = CondyleFit_i2 ;
        DiagInfo.rotationMatrices.CS3 = R_CS3 ;
        
    end
    %% Coordinate system Miranda et al 2010
    
    miranda_epicondylePoint_1 = computeIntersectionLineSurface(V,axisCondyles_i2,pointCondyles_i2' + (axisCondyles_i2.*20)',0);
    miranda_epicondylePoint_2 = computeIntersectionLineSurface(V,axisCondyles_i2,pointCondyles_i2' - (axisCondyles_i2.*20)',1);
    
    miranda_origin = mean([ miranda_epicondylePoint_1;miranda_epicondylePoint_2]) ; % centroid of the cylinder fit to the condyles
    miranda_ML = axisCondyles_i2; % Zdir
    miranda_AP = cross(miranda_ML,IA_Diaphysis(:,1)) ;
    miranda_AP = miranda_AP / rssq(miranda_AP) ;
    miranda_PD = cross(miranda_ML,miranda_AP) ;
    miranda_PD = miranda_PD / rssq(miranda_PD) ;
    strut_FrameM.miranda_origin = miranda_origin;
    strut_FrameM.miranda_ML = miranda_ML;
    strut_FrameM.miranda_AP = miranda_AP;
    strut_FrameM.miranda_PD = miranda_PD;
%     save('D:\project\4DTTTG\data\001\strut_FrameM', 'strut_FrameM'); 
    if nargout > 1
        DiagInfo.miranda.epicondyle_1 = miranda_epicondylePoint_1 ;
        DiagInfo.miranda.epicondyle_2 = miranda_epicondylePoint_2 ;
        
        mirandaCS.origin = miranda_origin;
        mirandaCS.X = miranda_AP';
        mirandaCS.Y = miranda_PD';
        mirandaCS.Z = miranda_ML';
        DiagInfo.miranda.CS = mirandaCS ;
    end
    
    %% Determining articulating surfaces
    %this is done by selecting vertices on the basis of 4 measures
    
    % we could evaluate the measures on different V's (e.g.
    % V_condyles2_CS4) which would simplify i.e. the plane distance metrics
    % we don't because of historical reasons
    % wich is to say I had this and it works so lets leave it
    
    %calculate the values of the 4 measures
    %angle between the miranda ML axis and the vertex normal
    vertNormDirection = abs(dot(repmat(miranda_ML', [size(V_Condyles_it2,1),1]),vertexnormals(F_Condyles_it2,V_Condyles_it2),2)) ;
    
    %distance to the condyle cutting plane
    backPlane = [pt1 axisCondyles_i2' , pt2-pt1] ;
    backPlaneNormal = cross(backPlane(:,4:6), backPlane(:, 7:9), 2) ;
    backPlaneNormal = backPlaneNormal / rssq(backPlaneNormal) ;
    distBackPlane = abs(sum(bsxfun(@times,backPlaneNormal , bsxfun(@minus, backPlane(:,1:3), V_Condyles_it2)), 2));
    
    %finding the plane seperating the condyles (this can be simplified i
    %think, why are we using CS4 here?)
    midPlaneAx1 = cross(miranda_ML,[0 0 1]')/rssq(cross(miranda_ML,[0 0 1]'));
    midPlaneAx2 = cross(miranda_ML,midPlaneAx1)/rssq(cross(miranda_ML,midPlaneAx1));
    R_SplitCondyles_CS4 = eye(3)*[midPlaneAx1 midPlaneAx2 miranda_ML]';
    % Rotate verticies in new CS (3rd axis miranda_MLaxis)
    miranda_origin_CS4 = (R_SplitCondyles_CS4*miranda_origin')';
    midPlane = reshape((R_SplitCondyles_CS4\[0 0 miranda_origin_CS4(3); 1 0 0; 0 1 0 ]'),[1,9]) ;
    
    midPlaneNormal = cross(midPlane(:,4:6), midPlane(:, 7:9), 2) ;
    midPlaneNormal = midPlaneNormal / rssq(midPlaneNormal) ;
    distMidPlane = abs(sum(bsxfun(@times,midPlaneNormal , bsxfun(@minus, midPlane(:,1:3), V_Condyles_it2)), 2));
    
    %threshold the metrics
    ccS(2,:)= (vertNormDirection < P_vertNormTol);
    ccS(3,:) = (distBackPlane > max(distBackPlane)*P_backDistFacTol)  ;
    ccS(4,:) = (distMidPlane > max(distMidPlane)*P_midDistFacTol)   ;
    
    %this does NOT take into account hte vertex normal measure because
    %   we didn't tune the paramters for that case.
    Ipost = all(ccS(3:4,:))' ;
    %distance to miranda's ML axis, normalised using points allowed by
    %previous measures (hence the "post")
    distAxPostPlanes = ones(size(distMidPlane)) ;
    distAxPostPlanes(Ipost) = bsxfun(@rdivide, vectorNorm3d(vectorCross3d(miranda_ML', bsxfun(@minus,miranda_origin, V_Condyles_it2(Ipost,:)))),vectorNorm3d(miranda_ML'));
    
    %threshold the last one
    ccS(5,:) = (distAxPostPlanes > max(distAxPostPlanes)*P_axDistPostFacTol)   ;
    
    %combine them
    ccS(1,:) = all(ccS(2:end,:)) ;
    %apply them
    V_artSurfaces = V_Condyles_it2(ccS(1,:)~=0,:) ;
    
    %seperate the two condyles
    sideTeller = sign(sum(bsxfun(@times,miranda_ML',bsxfun(@minus, V_artSurfaces , miranda_origin)),2)) ;
    side2Sign = sign(dot(miranda_ML', miranda_epicondylePoint_1 - miranda_origin)) ;
    V_artSurf_2 = V_artSurfaces(sideTeller==side2Sign,:) ;
    V_artSurf_1 = V_artSurfaces(sideTeller==-side2Sign,:) ;
    
    if nargout > 1
        DiagInfo.sections.artSurfaces.measureValues.midDist = distMidPlane ;
        DiagInfo.sections.artSurfaces.measureValues.vertexNormalAngle = vertNormDirection ;
        DiagInfo.sections.artSurfaces.measureValues.delimDist = distBackPlane ;
        DiagInfo.sections.artSurfaces.measureValues.distAx = distAxPostPlanes ;
        DiagInfo.sections.artSurfaces.backPlane = backPlane ;
        DiagInfo.sections.artSurfaces.midPlane =midPlane ;
        
        
        DiagInfo.sections.artSurfaces.constraintIndices = ccS ;
        DiagInfo.rotationMatrices.CS4 = R_SplitCondyles_CS4 ;
    end
    
    %% Fitting cylinder to the articulating points
    [artSurfCylAxisPoint_1, artSurfCylAxis_1, artSurfCylRadius_1, ...
        artSurfFitInfo_1.d, artSurfFitInfo_1.sigmah, artSurfFitInfo_1.conv, ...
        artSurfFitInfo_1.Vx0n, artSurfFitInfo_1.Van, artSurfFitInfo_1.urn, ...
        artSurfFitInfo_1.GNlog, artSurfFitInfo_1.a, artSurfFitInfo_1.R0, ...
        artSurfFitInfo_1.R] = lscylinder(V_artSurf_1, pointCondyles_i2, miranda_ML', radiusCondyles_i2, P_tolp(3), P_tolg(3));
    
    %find the cylinder center by projecting the surface CoM to the cylinder
    %axis
    VI_artSurfaces = (ccS(1,:)~=0) ;
    goodFi = all(reshape(any(bsxfun(@eq,F_Condyles_it2(:),find(VI_artSurfaces)),2),size(F_Condyles_it2)),2) ;
    
    Vi_artSurf_1 = VI_artSurfaces ;
    Vi_artSurf_1(Vi_artSurf_1==1) = sideTeller==-side2Sign ;
    Vofs = cumsum(~Vi_artSurf_1) ;
    
    F_artSurf_1 = F_Condyles_it2(goodFi,:) ;
    F_artSurf_1 = F_artSurf_1 - Vofs(F_artSurf_1) ;
    F_artSurf_1 = F_artSurf_1(all(F_artSurf_1>0,2),:) ;
    
    %compute the surface center of mass, and project to cylinder axis
    triCen = computeSurfaceCenterOfMass(F_artSurf_1,V_artSurf_1) ;
    answ = dot(triCen'-artSurfCylAxisPoint_1, artSurfCylAxis_1-artSurfCylAxisPoint_1) ;
    artSurfCenter_1 = artSurfCylAxisPoint_1 + artSurfCylAxis_1 / answ ;
    
    %the division above i dont get. but the assertion below tells me its right
    %     assert(rssq(cross((artSurfCenter_1' - triCen)/rssq(artSurfCenter_1' - triCen) ,artSurfCylAxis_1))- 1<eps(1)*100)
    if nargout > 1
        DiagInfo.fitCylinders.condyle1.center = artSurfCylAxisPoint_1 ;
        DiagInfo.fitCylinders.condyle1.axis = artSurfCylAxis_1 ;
        DiagInfo.fitCylinders.condyle1.radius = artSurfCylRadius_1 ;
        DiagInfo.fitCylinders.condyle1.fitInfo = artSurfFitInfo_1 ;
        
        DiagInfo.sections.articulatingSurface1.F = F_artSurf_1 ;
        DiagInfo.sections.articulatingSurface1.V = V_artSurf_1 ;
        DiagInfo.sections.articulatingSurface1.surfaceCentroid = triCen ;
        DiagInfo.sections.articulatingSurface1.centroid = artSurfCenter_1 ;
    end
    
    [artSurfCylAxisPoint_2, artSurfCylAxis_2, artSurfCylRadius_2, ...
        artSurfFitInfo_2.d, artSurfFitInfo_2.sigmah, artSurfFitInfo_2.conv, ...
        artSurfFitInfo_2.Vx0n, artSurfFitInfo_2.Van, artSurfFitInfo_2.urn, ...
        artSurfFitInfo_2.GNlog, artSurfFitInfo_2.a, artSurfFitInfo_2.R0, ...
        artSurfFitInfo_2.R] = lscylinder(V_artSurf_2, pointCondyles_i2, miranda_ML', radiusCondyles_i2, P_tolp(3), P_tolg(3));
    
    %same thing for the other articulating surface
    Vi_artSurf_2 = VI_artSurfaces ;
    Vi_artSurf_2(Vi_artSurf_2==1) = sideTeller==side2Sign ;
    Vofs = cumsum(~Vi_artSurf_2) ;
    
    F_artSurf_2 = F_Condyles_it2(goodFi,:) ;
    F_artSurf_2 = F_artSurf_2 - Vofs(F_artSurf_2) ;
    F_artSurf_2 = F_artSurf_2(all(F_artSurf_2>0,2),:) ;
    
    triCen = computeSurfaceCenterOfMass(F_artSurf_2,V_artSurf_2) ;
    answ = dot(triCen'-artSurfCylAxisPoint_2, artSurfCylAxis_2-artSurfCylAxisPoint_2) ;
    artSurfCenter_2 = artSurfCylAxisPoint_2 + artSurfCylAxis_2 / answ ;
    
    if nargout > 1
        DiagInfo.fitCylinders.condyle2.center = artSurfCylAxisPoint_2 ;
        DiagInfo.fitCylinders.condyle2.axis = artSurfCylAxis_2 ;
        DiagInfo.fitCylinders.condyle2.radius = artSurfCylRadius_2 ;
        DiagInfo.fitCylinders.condyle2.fitInfo = artSurfFitInfo_2 ;
        
        DiagInfo.sections.articulatingSurface2.F = F_artSurf_2 ;
        DiagInfo.sections.articulatingSurface2.V = V_artSurf_2 ;
        DiagInfo.sections.articulatingSurface2.surfaceCentroid = triCen ;
        DiagInfo.sections.articulatingSurface2.centroid = artSurfCenter_2 ;
    end
    
    %% Constructing final ML axis and origin using centerpoints of the cylinder fits
    
    % New ML axis
    MLaxis_new_GN = (artSurfCenter_2-artSurfCenter_1) ;%/(sqrt(sum((center_Lcylinder_GNax-center_Mcylinder_GNax).^2)));
    MLaxis_new_GN = MLaxis_new_GN /rssq(MLaxis_new_GN) ;
    origin_GN = mean([artSurfCenter_1,artSurfCenter_2]')  ; %#ok
    
    %% Create Lianne Coordinate system
    Lianne_originGN = origin_GN; % centroid of the cylinder fit to the condyles
    Lianne_MLaxisGN = MLaxis_new_GN; % Zdir
    Lianne_APaxisGN = cross(Lianne_MLaxisGN,IA_Diaphysis(:,1)) ;%Xdir
    Lianne_APaxisGN = Lianne_APaxisGN / rssq(Lianne_APaxisGN) ;
    Lianne_PDaxisGN = cross(Lianne_MLaxisGN,Lianne_APaxisGN); % Ydir
    Lianne_PDaxisGN = Lianne_PDaxisGN / rssq(Lianne_PDaxisGN) ;
    
    %% Apply the axis estimate to find the direciotn of the axes;
    femCoords.X= sign(dot(axEst.X,Lianne_APaxisGN)) * Lianne_APaxisGN';
    femCoords.Y = sign(dot(axEst.Y,Lianne_PDaxisGN)) * Lianne_PDaxisGN';
    femCoords.Z =sign(dot(axEst.Z,Lianne_MLaxisGN)) * Lianne_MLaxisGN';
    femCoords.origin = Lianne_originGN;
    
    %% Fix orientation of everything in fitCylinders to be readable when showing the struct in the console
    if nargout > 1
        DiagInfo.refFrame = femCoords ;
        
        DiagInfo.failed = 0 ;
        % fix the orientation of FitCylinder fields
        fcNs = fieldnames(DiagInfo.fitCylinders) ;
        for fii = 1:numel(fcNs)
            fn = fieldnames(DiagInfo.fitCylinders.(fcNs{fii}));
            for fiii = 1 : numel(fn) ;
                if all(size(DiagInfo.fitCylinders.(fcNs{fii}).(fn{fiii}))== [3,1])
                    DiagInfo.fitCylinders.(fcNs{fii}).(fn{fiii}) = DiagInfo.fitCylinders.(fcNs{fii}).(fn{fiii})';
                end
            end
        end
    end
catch emAll
    % For instance to short diaphysis causes the method to crash
    %this block makes sure the output still is somewhat sensible
    emAll.getReport
    
    DiagInfo.failed = 1;
    DiagInfo.errorMessage = emAll ;
    
    %also add phony data to mirandaCS
    if ~isfield(DiagInfo,'mirandaCS') ;
        mirandaCS.X =[1 0 0] ;
        mirandaCS.Y = [-1 0  0] ;
        mirandaCS.Z = [0 1 1] ;
        mirandaCS.origin = [ 0 0 0 ] ;
        DiagInfo.miranda.CS = mirandaCS ;
    end
    
    %phony data to reveal it goes wrong w/o crashing
    femCoords.X = [1 0 0] ;
    femCoords.Y = [-1 0  0] ;
    femCoords.Z = [0 1 1] ;
    femCoords.origin = [ 0 0 0 ] ;
    DiagInfo.refFrame = femCoords ;
end
end