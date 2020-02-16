function [ ] = drawSTL_TT_TG( FemurSTL, TibiaSTL, PatellaSTL, ...
    pointTT,pointTG,pointRL)

% draw the first time frame
%% Combine the bones together
BoneFace = [FemurSTL.F;PatellaSTL.F+size(FemurSTL.V,1);TibiaSTL.F+size(FemurSTL.V,1)+size(PatellaSTL.V,1)];
BoneVertices = [FemurSTL.V; PatellaSTL.V; TibiaSTL.V];

diffgap = 50;
xmax = max(BoneVertices(:,1))+diffgap;
xmin = min(BoneVertices(:,1))-diffgap;
ymax = max(BoneVertices(:,2))+diffgap;
ymin = min(BoneVertices(:,2))-diffgap;
zmax = max(BoneVertices(:,3))+diffgap;
zmin = min(BoneVertices(:,3))-diffgap;
figure;
%% Start
hold on;
axis equal;
axis off;%axis tight;
axis([xmin,xmax,ymin,ymax,zmin,zmax]);
hBone = patch('Faces', BoneFace, 'Vertices', BoneVertices, 'FaceColor', 1.0*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
light; labxyz;

set(hBone,'Vertices', BoneVertices);

scatter3(pointTT(1),pointTT(2),pointTT(3),'MarkerEdgeColor','R',  'MarkerFaceColor',[0 1 0]);
scatter3(pointTG(1),pointTG(2),pointTG(3),'MarkerEdgeColor','R',  'MarkerFaceColor','m');

scatter3(pointRL(1,1),pointRL(1,2),pointRL(1,3),'MarkerEdgeColor','B',  'MarkerFaceColor',[0 .75 .75]);
scatter3(pointRL(2,1),pointRL(2,2),pointRL(2,3),'MarkerEdgeColor','B',  'MarkerFaceColor',[0 .75 .75]);
plot3(pointRL(:,1),pointRL(:,2),pointRL(:,3),'k','LineWidth',2);

[ ProjPointTG ] = proj_point2line( pointRL(1,:),pointRL(2,:),pointTG );


[ ProjPointTT ] = proj_point2line( pointRL(1,:),pointRL(2,:),pointTT);
projPoints =[ProjPointTG;ProjPointTT];

plot3(projPoints(:,1),projPoints(:,2),projPoints(:,3),'color',[1 0 0],'LineWidth',10);
view(0,0);
end

