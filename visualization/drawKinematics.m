function [ ] = drawKinematics( FemurSTL, TibiaSTL, PatellaSTL, FemurRMs, TibiaRMs, PatellaRMs,...
    FemurCoordsPosDir,TibiaCoordsPosDir,PatellaCoordsPosDir,...
    pointTT,pointTG,pointTG_PM,pointTG_PB,pointRL,TT_Compare,TG_Flag)
%drawkinematics Summary of this function goes here
%   drawKinematics is aim to provide the kinematics anamation
%   FemurSTL, TibiaSTL, PatellaSTL:  the faces and vertices for femur,
%       tibia and patella  
%   FemurRMs, TibiaRMs, PatellaRMs:  the rigid rotation matrix including
%       the translation. The form is 4*4*N, which N is the amount of the
%       data. If N=10, then the dynamics data own 11 time frames. The
%       format is [U V W P]. So the coordinate transform is right multiply.
%   FemurCoordsPosDir,TibiaCoordsPosDir,PatellaCoordsPosDir: the origin
%       position and the direction of the coordinate 
%   pointTT: the coordinate of tibial T
%   pointTG: the coordinate of trogh groove
%   pointRL: the reference line on posterior condyler, Row 1 is closer to
%   the medial direction, while Row 2 is closer to lateral direction

%   Hao 2017.08-19

% draw the first time frame
%% Combine the bones together
BoneFace = [FemurSTL.faces;PatellaSTL.faces+size(FemurSTL.vertices,1);TibiaSTL.faces+size(FemurSTL.vertices,1)+size(PatellaSTL.vertices,1)];
BoneVertices = [FemurSTL.vertices; PatellaSTL.vertices; TibiaSTL.vertices];
diffgap = 50;
xmax = max(BoneVertices(:,1))+diffgap;
xmin = min(BoneVertices(:,1))-diffgap;
ymax = max(BoneVertices(:,2))+diffgap;
ymin = min(BoneVertices(:,2))-diffgap;
zmax = max(BoneVertices(:,3))+diffgap;
zmin = min(BoneVertices(:,3))-diffgap;
struct_MinMax.xmin = xmin;
struct_MinMax.xmax = xmax;
struct_MinMax.ymin = ymin;
struct_MinMax.ymax = ymax;
struct_MinMax.zmin = zmin;
struct_MinMax.zmax = zmax;

%% Start
nTF_Num=size(FemurRMs,3);
num = 1;
nstep =1;
figure;
TTTGdis = zeros(11,1);
while(num<=nTF_Num)
    if isempty(get(0,'children'))
        return
    end
    cla;
    axis equal;
    axis off;%axis tight;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]);
    hBone = patch('Faces', BoneFace, 'Vertices', BoneVertices, 'FaceColor', 1.0*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    light; labxyz;
    TFormFem = affine3d(FemurRMs(:,:,num));
    TFormPat = affine3d(PatellaRMs(:,:,num));
    TFormTib = affine3d(TibiaRMs(:,:,num));
    BoneVertices = [FemurSTL.vertices;transformPointsInverse(TFormFem,transformPointsForward(TFormPat, PatellaSTL.vertices));transformPointsInverse(TFormFem,transformPointsForward(TFormTib, TibiaSTL.vertices))];
    set(hBone,'Vertices', BoneVertices);
%     draw_mesh_coord(BoneFace,BoneVertices,FemurCoordsPosDir,TibiaCoordsPosDir);
    
    hp = plotCoordsNew(FemurCoordsPosDir(1:3, 1:3),  FemurCoordsPosDir(1:3, 4)',struct_MinMax);
    hp = plotCoordsNew((PatellaCoordsPosDir(1:3, 1:3)'*PatellaRMs(1:3,1:3,num)*FemurRMs(1:3, 1:3,num)')', transformPointsInverse(TFormFem,transformPointsForward(TFormPat, PatellaCoordsPosDir(1:3, 4)')),struct_MinMax);
    hp = plotCoordsNew((TibiaCoordsPosDir(1:3, 1:3)'*TibiaRMs(1:3,1:3,num)*FemurRMs(1:3, 1:3,num)')', transformPointsInverse(TFormFem,transformPointsForward(TFormTib, TibiaCoordsPosDir(1:3, 4)')),struct_MinMax);
%     tP = text(hBone(1).XData(1), hBone(1).YData(1), hBone(1).ZData(1), sprintf('P%d', num), 'FontSize', 12);
    pointTT_T = transformPointsInverse(TFormFem,transformPointsForward(TFormTib, pointTT));
    pointTG_PM_T = transformPointsInverse(TFormFem,transformPointsForward(TFormPat, pointTG_PM));
    pointTG_PB_T = transformPointsInverse(TFormFem,transformPointsForward(TFormPat, pointTG_PB));
    if (~isempty(TT_Compare))    
        pointTT_TNew = transformPointsInverse(TFormFem,transformPointsForward(TFormTib, TT_Compare));
        scatter3(pointTT_TNew(1),pointTT_TNew(2),pointTT_TNew(3),'MarkerEdgeColor','R',  'MarkerFaceColor',[1 0 0]);
    end
    scatter3(pointTT_T(1),pointTT_T(2),pointTT_T(3),'MarkerEdgeColor','k',  'MarkerFaceColor',[0 0 1]);
    scatter3(pointTG(1),pointTG(2),pointTG(3),'MarkerEdgeColor','r',  'MarkerFaceColor',[1 0 0]);
%     scatter3(pointTG_PM_T(1),pointTG_PM_T(2),pointTG_PM_T(3),'MarkerEdgeColor','R',  'MarkerFaceColor',[0 1 0]);
%     scatter3(pointTG_PB_T(1),pointTG_PB_T(2),pointTG_PB_T(3),'MarkerEdgeColor','R',  'MarkerFaceColor',[0 0 1]);
    scatter3(pointRL(1,1),pointRL(1,2),pointRL(1,3),'MarkerEdgeColor','B',  'MarkerFaceColor',[0 .75 .75]);
    scatter3(pointRL(2,1),pointRL(2,2),pointRL(2,3),'MarkerEdgeColor','B',  'MarkerFaceColor',[0 .75 .75]);
    plot3(pointRL(:,1),pointRL(:,2),pointRL(:,3));
    switch TG_Flag
        case 1
            [ ProjPointTG ] = proj_point2line( pointRL(1,:),pointRL(2,:),pointTG );
        case 2
            [ ProjPointTG ] = proj_point2line( pointRL(1,:),pointRL(2,:),pointTG_PM_T );
        case 3
            [ ProjPointTG ] = proj_point2line( pointRL(1,:),pointRL(2,:),pointTG_PB_T );
    end
    
    [ ProjPointTT ] = proj_point2line( pointRL(1,:),pointRL(2,:),pointTT_T );
    projPoints =[ProjPointTG;ProjPointTT];
    directionRL = projPoints(2,:)-projPoints(1,:);
    signFlag =sign(dot(FemurCoordsPosDir(1:3,3),directionRL));
    if(signFlag<0)
        plot3(projPoints(:,1),projPoints(:,2),projPoints(:,3),'color',[0 0 0],'LineWidth',10);
%         TTTGdis(num+1) = -sum((ProjPointTG-ProjPointTT).^2).^0.5;
    else
        plot3(projPoints(:,1),projPoints(:,2),projPoints(:,3),'color',[0 0 0],'LineWidth',10);
%         TTTGdis(num+1) = sum((ProjPointTG-ProjPointTT).^2).^0.5;
    end
    
    drawnow;
    pause(0.1);
    legend('off')
%     view(0,0)
%     str_file = ['F:\work\ppt\DKCT\figure\dynamic_coronal' num2str(num) '.png'];
%     set(gcf,'color','w');
%     print(str_file,'-dpng','-r300')
%     export_fig(gcf,str_file, '-r300');
    num=num+nstep;    
    if(num>nTF_Num)
        num=1;
    end
end

end

