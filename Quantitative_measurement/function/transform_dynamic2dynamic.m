function [DynamicCoordinateNew] = transform_dynamic2dynamic(OriginalCoordinate,mat_affine3d)
%transform_dynamic2dynamic is aimed to transform the 3d point from dynamic
%image based on the rotation matrix of CPD
%   Detailed explanation goes here
% Input:
% OriginalCoordinate: the original coordinate with form:
%         u 0     Xu Yu Zu 0
%         v 0  =  Xv Yv Zv 0
%         w 0     Xw Yw Zw 0 
%         o 1     Xo Yo Zo 1
% mat_affine3d: the rotation matrix
%
% Output:
% DynamicCoordinateNew: the rotated coordinate with the same form of OriginalCoordinate

% Hao
% 2018-09-16

DynamicCoordinateNew = OriginalCoordinate;
DynamicCoordinateNew(4,1:3) = transformPointsForward(mat_affine3d,DynamicCoordinateNew(4,1:3));
rotate = mat_affine3d;
rotate.T(4,1:3) = 0; 
DynamicCoordinateNew(1:3,1:3) = transformPointsForward(rotate,OriginalCoordinate(1:3,1:3));
end

