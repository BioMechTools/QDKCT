function  plotCylinder(MpointAxis,LpointAxis,dirAxis,radius,ccCyl)
%GENERATEPLOTABLECYLINDERPOINTS generates inputs to surf for visualising a
%cylinder
% 
% Input:
% MpointAxis: medial point cylinderaxis
% LpointAxis: lateral point cylinderaxis
% dirAxis: direction vector cylinder
% radius: cylinder radius
%
% Output:
% plot_cylinder
%
% Example use of output:
%
% surf(plot_cylinder(:,:,1),plot_cylinder(:,:,2),plot_cylinder(:,:,3),...
%    'FaceColor','c','EdgeColor','none','FaceAlpha',0.5,...
%    'SpecularColorReflectance', 0, 'SpecularExponent', 50,'DiffuseStrength',1);
%
% ORL Nijmegen, 2012

%robustify against odd input orientation
MpointAxis = reshape(MpointAxis,[1 3]) ;
LpointAxis = reshape(LpointAxis,[1 3]) ;
dirAxis = reshape(dirAxis,[3 1]) ;

% Rotation matrix 
R(:,3) = dirAxis;
R(:,1) = cross(dirAxis,cross(dirAxis,[0;1;0])); R(:,1)  = R(:,1) /rssq(R(:,1) ) ;
R(:,2) = cross(dirAxis,R(:,1)); R(:,2) = R(:,2) / rssq(R(:,2)) ;
rot = eye(3)*R';

pointAxis1R = rot*MpointAxis';
pointAxis2R = rot*LpointAxis';

% Fit circle
[plot_circleX1R,plot_circleY1R] = plotcircle2D_nested(pointAxis1R(1),pointAxis1R(2),radius);
plot_circle1R = [plot_circleX1R;plot_circleY1R;zeros(size(plot_circleX1R))+pointAxis1R(3)];

[plot_circleX2R,plot_circleY2R] = plotcircle2D_nested(pointAxis2R(1),pointAxis2R(2),radius);
plot_circle2R = [plot_circleX2R;plot_circleY2R;zeros(size(plot_circleX2R))+pointAxis2R(3)];

% Project back
plot_circle1 = (inv(rot)*plot_circle1R);  % points on cylinder
plot_circle2 = (inv(rot)*plot_circle2R);

plot_cylinder(1,:,:) = plot_circle1';
plot_cylinder(2,:,:) = plot_circle2';

    surf(plot_cylinder(:,:,1),plot_cylinder(:,:,2),plot_cylinder(:,:,3),...
        'FaceColor',ccCyl,'FaceAlpha',0.5,'EdgeAlpha',0,...
        'SpecularColorReflectance', 0, 'SpecularExponent', 50,'DiffuseStrength',1);

function [xp, yp] = plotcircle2D_nested(xc,yc,radius)
%PLOTCIRCLE2D x and y position of points on a circle
% xc and yc are the coordinates of the center of the circle
% r is the radius of the circle
ang = 0:0.05:2*pi+0.05;
xp = xc + radius*cos(ang);
xp = xp(:)';
yp = yc + radius*sin(ang);
yp = yp(:)';