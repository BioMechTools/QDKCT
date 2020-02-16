function [dis_vertices,dist_surf,dist_vertices] = registration_error(T_V,T_F, O_V,O_F, C_list)
%% this function is aimed to caluclate the point to surface distance based on restered object and target object
%%% 
% Input:
%   T_V: vertices of target
%   T_F: faces of target
%   O_V: vertices of registered object
%   O_F: faces of registered object
%   C_list: correspondece point from target to register 
diff = O_V - T_V(C_list,:);
dis_vertices = sqrt(sum(diff.^2,2));
TR = triangulation(O_F,O_V);
% VN = vertexNormal(TR);
FV.faces = T_F;
FV.vertices = T_V;
[dist_surf,dist_vertices] = point2trimesh(FV, 'QueryPoints', O_V); 
% [dist_surf,~,~,~,~]=raysurf(O_V,VN,,);
end