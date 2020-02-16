function [vertices_cut,faces_cut] = get_cut_bone_static(vertices_origin,faces_origin, vertices_target,faces_target,str_bone)
%get_cut_bone Summary of this function goes here
%Input:
%   F_V,F_F: vertices and faces of femur
%   T_V,T_F: vertices and faces of tibia
%   P_V,P_F: vertices and faces of patella
%   str_bone: 'femur' or tibia, bone type that will be cut
%
%   Hao
%   2018-10-01


IA = pca(vertices_origin) ; % calculate the principle axis for the vertices
v_origin_proj = vertices_origin *IA(:,1);
len_origin = max(v_origin_proj(:))-min(v_origin_proj(:));
%%% length of the target bone
V_proj_target = vertices_target *IA(:,1);
len_target = max(V_proj_target(:))-min(V_proj_target(:));
switch str_bone
    case 'femur'
        percRange = [0 round(len_target/len_origin*100)];
    case 'tibia'
        percRange = [100-round(len_target/len_origin*100),100];
end
[refL,VcsL,IAtemp] = computeReferenceLength(vertices_origin) ;
[F_tnew, V_t_new1] = cutMeshByPercent(faces_origin,VcsL,percRange,1,1) ;

V_t_new1 = (IAtemp'\V_t_new1')' ;
figure;plotsurf(V_t_new1,F_tnew);
faces_cut = F_tnew;
vertices_cut = V_t_new1;

end

