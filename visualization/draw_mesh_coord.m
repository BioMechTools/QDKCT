function [] = draw_mesh_coord(BoneFace,BoneVertices,femCoords,tibCoords)
%draw_mesh_coord Summary of this function goes here
%   Detailed explanation goes here
figure;
axis equal;axis off;
hBone = patch('Faces', BoneFace, 'Vertices', BoneVertices, 'FaceColor', 1.0*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
light; labxyz;
set(hBone,'Vertices', BoneVertices);
% FemurCoordsPosDir1 = [femCoords.AP femCoords.PD femCoords.ML];
% hp = plotCoords(FemurCoordsPosDir1, femCoords.Origin);
% TibiaCoordsPosDir = [tibCoords.AP tibCoords.PD tibCoords.ML];
% hp = plotCoords(TibiaCoordsPosDir, (tibCoords.Origin));
drawnow;
end

