function [Epsilon,theta,D,d] = get_para_Yao(staticLandmarks,CT_Frame,tibCoords)
%get_para_Yao Summary of this function goes here
%   Detailed explanation goes here

%%% unit direction for PCL
dir_PCL = (staticLandmarks(2,:) - staticLandmarks(1,:))/norm(staticLandmarks(2,:) - staticLandmarks(1,:));
% if str_side == 'L'
dir_ML_CT = CT_Frame.ML;%[1 0 0];
% else
%     dir_ML_CT = -CT_Frame.ML;
% end
[Epsilon] = get_angle_vectors(dir_PCL,dir_ML_CT);
%%% D d %% TT is the original point
vec_TT_TG = staticLandmarks(3,:) - staticLandmarks(4,:);
[angle] = get_angle_vectors(vec_TT_TG,CT_Frame.PD);%[0 0 1]
D = norm(vec_TT_TG)*cosd(angle);
[angle] = get_angle_vectors(vec_TT_TG,CT_Frame.AP);%[0 1 0]
d = norm(vec_TT_TG)*cosd(angle);
% PD = [0 0 1];
% ML = [1 0 0];
A = [CT_Frame.PD' CT_Frame.ML'];
[v] = get_proj_vector_plane(tibCoords.PD',A);
[theta] = get_angle_vectors(v',CT_Frame.PD);
if v(1)<0
    theta = -theta;
end
end

function [angle_cos] = get_angle_vectors(u,v)
CosTheta = dot(u,v)/(norm(u)*norm(v));
angle_cos = acosd(CosTheta);
end

function [v] = get_proj_vector_plane(u,A)
%%%
% A is 3*2 matrix, each column is a direction of the plane
% u is 3*1 vector, and is the original vector
% v is 3*1 vector, and is the projected vector
v = A*inv(A'*A)*A'*u;

end