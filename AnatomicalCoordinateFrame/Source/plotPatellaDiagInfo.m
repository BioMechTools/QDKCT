function plotPatellaDiagInfo( DiagInfo,stlData )
%PLOTPATELLADIAGINFO Visualised the different steps to undertaken to
%generate the patellar reference frame in CSpatellaMAX
%   plotPatellaDiagInfo(DiagInfo, F, V)
%
%---Input
%   DiagInfo    -   Second output of CSpatellaMAX. 
%   F           -   Faces. N by 3 integer array
%   V           -   Vertices. N by 3 array.
%       Make sure DiagInfo is generated on the basis of these F and V
%
%  Example usage:
%    [PatellaFrame, DiagInfo] = CSpatellaMAX(F, V, XYZest) ;
%    plotPatellaDiagInfo(DiagInfo,F,V) ;
%
%  (V1.1) ORL Nijmegen, Max Bakker 2016
%
%  see also: CSpatellaMAX 

%Chagelog
%1.1, fixed a typo with the getReport at the end

%try block makes sure no error is thrown when the ERCrefFramePatella didnt finish
F = stlData.F ;
V = stlData.V ;
try
    mncVec = mean([DiagInfo.refFrame.X ;DiagInfo.refFrame.Y;DiagInfo.refFrame.Z]) ;
    az0 =  tan(mncVec(2)/mncVec(1)) ;
    el0 = cos(rssq(mncVec([1,2])) / rssq(mncVec)) ;
    
    
    patellaCS = DiagInfo.refFrame ;
    
    figure(123) ; clf
    
    subplot(2,2,1)
    axis equal; hold on ;
    title('Patellar IA and Ant surf')
    patch('Faces',F,'Vertices',V,'FaceColor',rand(3,1)*.4,'FaceAlpha',0.4,'EdgeColor','none')
    patch('Faces',DiagInfo.sections.AnteriorPatella.F,'Vertices',DiagInfo.sections.AnteriorPatella.V,'FaceColor','none','EdgeAlpha',0.4)
    plotCoords(DiagInfo.IA.IA,DiagInfo.IA.CM,'Color',eye(3)*0.4)
    % ezplotm([DiagInfo.polePosition1.guess;DiagInfo.polePosition2.guess],'-k')
    view(az0,el0)
    
    
    subplot(2,2,2)
    axis equal; hold on ;
    title('fit plane to determine AP')
    patch('Faces',DiagInfo.sections.AnteriorPatella.F,'Vertices',DiagInfo.sections.AnteriorPatella.V,'FaceAlpha',0.5)
    hAp = drawPlane3d([DiagInfo.sections.AnteriorPatella.planeFitCentroid' patellaCS.Y patellaCS.Z]) ;
    hlAP = drawLine3d([DiagInfo.sections.AnteriorPatella.planeFitCentroid' patellaCS.X]) ;
    set(hlAP,'Color','r','LineWidth',2)
    set(hAp,'FaceAlpha',0.4)
    view(az0,el0)
    
    
    subplot(2,2,3)
    axis equal; hold on ;
    patch('Faces',F,'Vertices',V,'FaceColor','none','EdgeAlpha',0.5);
    title('Pole Finding');
    ccIPF = cool(numel(DiagInfo.polePosition)*2+4);
    ccIPF = ccIPF(4:2:end,:) ;
    hcm = ezplotm(patellaCS.origin,'k.','MarkerSize',20) ;
    for gi = 1 : numel(DiagInfo.polePosition) ;
        thisFitS = DiagInfo.polePosition(gi) ;
        hg =ezplotm(thisFitS.guess,'.') ;
        hl = drawLine3d([patellaCS.origin,patellaCS.origin-thisFitS.guess]) ;
        set(hl,'LineStyle',':','Color',ccIPF(gi,:),'LineWidth',2)
        set(hg,'MarkerSize',25,'Color',ccIPF(gi,:));
        view(az0,el0)
    end
    set(hg,'Color',[0 1 0]);
    set(hl,'Color',[0 1 0],'LineStyle','-');
    
    
    
    subplot(2,2,4)
    axis equal; hold on ;
    title('Patlla CS')
    patch('Faces',F,'Vertices',V,'FaceColor',rand(3,1)*.4,'FaceAlpha',0.4,'EdgeAlpha',0.3)
    plotCoords(patellaCS) ;
    thisFitS =  DiagInfo.polePosition(end) ;
    hl = drawLine3d([patellaCS.origin,patellaCS.origin-thisFitS.guess]) ;
    set(hl,'LineStyle','-','Color',[0 1  0],'LineWidth',2)
    ezplotm(patellaCS.origin,'k.','MarkerSize',14)
    hAp = drawPlane3d([DiagInfo.sections.AnteriorPatella.planeFitCentroid' patellaCS.Y patellaCS.Z]) ;
    set(hAp,'FaceAlpha',0.4)
    hlAP = drawLine3d([DiagInfo.sections.AnteriorPatella.planeFitCentroid' patellaCS.X]) ;
    set(hlAP,'Color','r','LineWidth',2)
    view(az0,el0)
    
    
    figure(1111) ;clf
    plotPlacer( numel(DiagInfo.polePosition) )
    scaleF = DiagInfo.param.infPoleFindScaleF(1) ;
    dotSz = 25 ;
    for gi = 1 : numel(DiagInfo.polePosition) ;
        scaleA = [1 scaleF 1] ;
        thisFitS = DiagInfo.polePosition(gi) ;
        V_find_CSl_scaled = bsxfun(@times,thisFitS.V_CSl,scaleA) ;
        
        plotPlacer ;
        axis equal; hold on ;
        Vscale_CS0 = (thisFitS.inAxes'\V_find_CSl_scaled')';
        cc = parula(64) ;
        cc  = cc(end:-1:1,:) ;
        scatter3(Vscale_CS0(:,1),Vscale_CS0(:,2),Vscale_CS0(:,3),abs(Vscale_CS0(:,3))*1+dotSz,thisFitS.distA,'.')
        colormap(cc) ;
        
        [~,maxi] = max(thisFitS.distA) ;
        title(['Pole Finding it:' num2str(gi) ' sf:' num2str(scaleF)]) ;
        
        view(az0,el0)
        labxyz
        ezplotm(Vscale_CS0(maxi,:),'mo','MarkerSize',20)
        ezplotm(Vscale_CS0(maxi,:),'mx','MarkerSize',20)
        ezplotm(zeros(1,3),'ko','MarkerSize',16)
    end
    
catch em
    if ~DiagInfo.failed
        rethrow(em) ;
    else
        DiagInfo.errorMessage.getReport ;
    end
end
end