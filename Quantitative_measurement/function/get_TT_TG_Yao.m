function [TT_TG_corrected] = get_TT_TG_Yao(Epsilon,theta,D,d,TT_TG_dist,str_side)
%get_para_Yao Summary of this function goes here
%   Detailed explanation goes here
% Epsilon = 180-Epsilon;
% theta = -theta;

% if Epsilon >90
%     Epsilon = -(180-Epsilon);
% end
% if theta >90
%     theta = -(180-theta);
% end
% if str_side == 'R'
%     Epsilon = -Epsilon;
%     theta = -theta;
% end
TT_TG_corrected = d*sind(Epsilon)*(1-cosd(theta)) + TT_TG_dist*cosd(theta) + D*cosd(Epsilon)*sind(theta);

end

