function [] = draw_TT(V_Tibia,F_Tibia,femCoords,tibCoords,TTPoint3D,str_figure,n_subj,str_side)
%draw_TT Summary of this function goes here
%   Detailed explanation goes here

draw_mesh_coord(F_Tibia,V_Tibia,femCoords,tibCoords);
axis tight
hold on;
if size(TTPoint3D,1)==1
    scatter3(TTPoint3D(1),TTPoint3D(2),TTPoint3D(3), 'MarkerEdgeColor','M','MarkerFaceColor',[0 .75 .75]);
else
    scatter3(TTPoint3D(1,1),TTPoint3D(1,2),TTPoint3D(1,3), 'MarkerEdgeColor','M','MarkerFaceColor','g');
    scatter3(TTPoint3D(2:4,1),TTPoint3D(2:4,2),TTPoint3D(2:4,3), 'MarkerEdgeColor','M','MarkerFaceColor','r');
end
hold off;
view(-9,7);
str_file = [str_figure num2str(n_subj) 'TT_' str_side '.png'];
saveas(gcf,str_file);
end

