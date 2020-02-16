function [intpt_axis] = computeIntersectionLineSurface(Vertices,dir_axis,pt_axis,int01)
%INTERSECTIONAXISSURFACE find the point of intersection between a line
% and a surface represented by vertices
%
% [intpt_axis] = computeIntersectionPointLineSurface(Vertices,dir_axis,pt_axis,int01)
%
% Input:
% Vertices = vertices of a surface
% dir_axis = direction cosine of the intersection axis
% pt_axis = point on the intersection axis
% int01 = 0 - intersection axis/surface in the direction of the line
%         1 - intersection axis/surface oposite to the direction of the line
%
% Output:
% intpt_axis = intersection point of the axis with the surface of the vertices
%
% ORL Nijmegen, WJ Zevenbergen June 2012
% Revised by Max November 2016

% Rotation matrix dir_axis
z_axis = [0 0 1]';          % z-axis
B = cross(dir_axis,z_axis);
if rssq(B)==0
    B = cross(dir_axis,[0 0.5 0.5]');
end
uv_B = B/rssq(B);
C = cross(dir_axis,uv_B);
uv_C = C/rssq(C);
R = [uv_B,uv_C,dir_axis]';  % applied rotation - 3rd rotation is dir-axis


% Rotation
pt_axis_r = R*pt_axis';     % rotation point on axis
V_r = (R*Vertices')';       % rotated vertices


V_r_plus = (V_r(:,3) > pt_axis_r(3)); % logicals -> position along the axis of the rotated vertices > point on the axis

% Compute the distance between the rotated vertices and point on the axis
% projected perpendicular to the axis direction
distInPlane = rssq(bsxfun(@minus,V_r(:,[1,2]),pt_axis_r([1,2])'),2)' ;
[~, ind_distance] = sort(distInPlane,'ascend') ;
[~,i] = min(V_r_plus(ind_distance)==int01) ;
intpt_axis_guess = Vertices(ind_distance(i),:) ;


% in stead of sorting the distance and finding hte first on the right side
% i here tried to find the minimum distance on the right side, and then
% correct the index offset. couldn't get it to work. 
% correctSideVertI = V_r_plus~=int01 ;
% [~,minI] = min(distInPlane(correctSideVertI)) ;
% % correct for correctSideVertI not having all vertices
% i2 = minI + sum(~correctSideVertI(1:minI)) ;
% intpt_axis_guess2 = Vertices(minI,:)

vect_ptaxis_intptaxis = intpt_axis_guess - pt_axis;
intpt_axis = pt_axis + (dir_axis.*dot(dir_axis,vect_ptaxis_intptaxis) )';
