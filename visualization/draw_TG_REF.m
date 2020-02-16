function [] = draw_TG_REF(V_Femur,F_Femur,femCoords,tibCoords,TGPoint3D,LinePoint1_3d,LinePoint2_3d,str_figure,n_subj,str_side)
%draw_TT Summary of this function goes here
%   Detailed explanation goes here

% draw_mesh_coord(F_Femur,V_Femur,femCoords,tibCoords);
figure;
axis equal;axis off;
hBone = patch('Faces', F_Femur, 'Vertices', V_Femur, 'FaceColor', 1.0*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
light; labxyz;
set(hBone,'Vertices', V_Femur);
axis tight
hold on;
if size(TGPoint3D,1)==1
    scatter3(TGPoint3D(1),TGPoint3D(2),TGPoint3D(3), 'MarkerEdgeColor','M','MarkerFaceColor',[0 .75 .75]);
    scatter3(LinePoint1_3d(1),LinePoint1_3d(2),LinePoint1_3d(3), 'MarkerEdgeColor','M','MarkerFaceColor',[0 .75 .75]);
    scatter3(LinePoint2_3d(1),LinePoint2_3d(2),LinePoint2_3d(3), 'MarkerEdgeColor','M','MarkerFaceColor',[0 .75 .75]);
elseif size(TGPoint3D,1)==4
    scatter3(TGPoint3D(1,1),TGPoint3D(1,2),TGPoint3D(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','g');
    scatter3(TGPoint3D(2:4,1),TGPoint3D(2:4,2),TGPoint3D(2:4,3), 'MarkerEdgeColor','M','MarkerFaceColor','r');
    
    scatter3(LinePoint1_3d(1,1),LinePoint1_3d(1,2),LinePoint1_3d(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','y');
    scatter3(LinePoint1_3d(2:4,1),LinePoint1_3d(2:4,2),LinePoint1_3d(2:4,3), 'MarkerEdgeColor','M','MarkerFaceColor','b');
    
    scatter3(LinePoint2_3d(1,1),LinePoint2_3d(1,2),LinePoint2_3d(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','y');
    scatter3(LinePoint2_3d(2:4,1),LinePoint2_3d(2:4,2),LinePoint2_3d(2:4,3), 'MarkerEdgeColor','M','MarkerFaceColor','b');
    
    line0 = [LinePoint1_3d(1,:);LinePoint2_3d(1,:)];
    line1 = [LinePoint1_3d(2,:);LinePoint2_3d(2,:)];
    line2 = [LinePoint1_3d(3,:);LinePoint2_3d(3,:)];
    line3 = [LinePoint1_3d(4,:);LinePoint2_3d(4,:)];
    plot3(line0(:,1),line0(:,2),line0(:,3),'LineWidth',1.5);
    plot3(line1(:,1),line1(:,2),line1(:,3),'LineWidth',1.5);
    plot3(line2(:,1),line2(:,2),line2(:,3),'LineWidth',1.5);
    plot3(line3(:,1),line3(:,2),line3(:,3),'LineWidth',1.5);
else
    scatter3(TGPoint3D(1,1),TGPoint3D(1,2),TGPoint3D(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','k');
    scatter3(TGPoint3D(2,1),TGPoint3D(2,2),TGPoint3D(2,3), 'MarkerEdgeColor','M','MarkerFaceColor','g');
    scatter3(TGPoint3D(3:5,1),TGPoint3D(3:5,2),TGPoint3D(3:5,3), 'MarkerEdgeColor','M','MarkerFaceColor','r');
    
    scatter3(LinePoint1_3d(1,1),LinePoint1_3d(1,2),LinePoint1_3d(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','m');
    scatter3(LinePoint1_3d(2,1),LinePoint1_3d(2,2),LinePoint1_3d(2,3), 'MarkerEdgeColor','M','MarkerFaceColor','y');
    scatter3(LinePoint1_3d(3:5,1),LinePoint1_3d(3:5,2),LinePoint1_3d(3:5,3), 'MarkerEdgeColor','M','MarkerFaceColor','b');
    
    scatter3(LinePoint2_3d(1,1),LinePoint2_3d(1,2),LinePoint2_3d(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','m');
    scatter3(LinePoint2_3d(2,1),LinePoint2_3d(2,2),LinePoint2_3d(2,3), 'MarkerEdgeColor','M','MarkerFaceColor','y');    
    scatter3(LinePoint2_3d(3:5,1),LinePoint2_3d(2:4,2),LinePoint2_3d(3:5,3), 'MarkerEdgeColor','M','MarkerFaceColor','b');
    
    line0 = [LinePoint1_3d(1,:);LinePoint2_3d(1,:)];
    line1 = [LinePoint1_3d(2,:);LinePoint2_3d(2,:)];
    line2 = [LinePoint1_3d(3,:);LinePoint2_3d(3,:)];
    line3 = [LinePoint1_3d(4,:);LinePoint2_3d(4,:)];
    line4 = [LinePoint1_3d(5,:);LinePoint2_3d(5,:)];
    plot3(line0(:,1),line0(:,2),line0(:,3),'LineWidth',1.5);
    plot3(line1(:,1),line1(:,2),line1(:,3),'LineWidth',1.5);
    plot3(line2(:,1),line2(:,2),line2(:,3),'LineWidth',1.5);
    plot3(line3(:,1),line3(:,2),line3(:,3),'LineWidth',1.5);
    plot3(line4(:,1),line4(:,2),line4(:,3),'LineWidth',1.5);
end
% scatter3(TGPoint3D(1),TGPoint3D(2),TGPoint3D(3), 'MarkerEdgeColor','M','MarkerFaceColor',[0 .75 .75]);
view(-9,7);
hold off;
str_file = [str_figure num2str(n_subj) 'TG_' str_side '.png'];
saveas(gcf,str_file);
end