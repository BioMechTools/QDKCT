function [TTTGdis] = cal_TT_TG_dis(pointRL,TG,TT,v_AP,str_BoneSide)
%cal_TT_TG_dis Summary of this function goes here
%   Detailed explanation goes here

[ ProjPointTG ] = proj_point2line( pointRL(1,:),pointRL(2,:),TG );
[ ProjPointTT ] = proj_point2line( pointRL(1,:),pointRL(2,:),TT);
%%% determine the sign positive direction to the project TG or not
projPoints =[ProjPointTG;ProjPointTT];
directionRL = projPoints(2,:)-projPoints(1,:);
signFlag =sign(dot(v_AP(1),directionRL(1)));
TTTGdis = sum((ProjPointTT-ProjPointTG).^2).^0.5;
if strcmp('R',str_BoneSide)
    if(signFlag<0)
        TTTGdis = -TTTGdis;
    end
else
    if(signFlag>0)
        TTTGdis = -TTTGdis;
    end
end
end

