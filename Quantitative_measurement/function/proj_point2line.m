function [ ProjPoint ] = proj_point2line( LinePA,LinePB,PointP )
%proj_point2line is aimed to project a point to a line
%   Detailed explanation goes here
ap = PointP-LinePA;
ab = LinePB-LinePA;
ProjPoint = LinePA + dot(ap,ab)/dot(ab,ab).* ab;
end

