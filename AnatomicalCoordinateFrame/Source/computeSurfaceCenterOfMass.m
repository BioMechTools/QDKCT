function [Centroid] = computeSurfaceCenterOfMass(F,V)
%computeSurfaceCenterOfMass.M centroid of a surface object
%
% Input:
% F - faces: contains the vertex lists defining each triangle face [n x 3]
% V - vertices: contains the vertices for all triangles [3*n x 3]
%
% Output:
% Centroid - centroid of a surface object
%
% ORL Nijmegen, WJ Zevenbergen May 2012

%% Step 1: compute surface area faces
Aface = zeros(size(F,1),1);
for i = 1:size(F,1)
    AB = [V(F(i,2),1)-V(F(i,1),1), V(F(i,2),2)-V(F(i,1),2), V(F(i,2),3)-V(F(i,1),3)];
    AC = [V(F(i,3),1)-V(F(i,1),1), V(F(i,3),2)-V(F(i,1),2), V(F(i,3),3)-V(F(i,1),3)];
    ABC = cross(AB,AC);
    Aface(i) = 0.5*(sqrt(ABC(1).^2+ABC(2).^2+ABC(3).^2));
end

%% Step 2: compute center of mass faces
CoMxx = zeros(1,size(F,1));
CoMyy = zeros(1,size(F,1));
CoMzz = zeros(1,size(F,1));
for i = 1:size(F,1)
    CoMxx(i) = (V(F(i,1),1)+V(F(i,2),1)+V(F(i,3),1))/3; % center of mass x-position
    CoMyy(i) = (V(F(i,1),2)+V(F(i,2),2)+V(F(i,3),2))/3; % center of mass y-position
    CoMzz(i) = (V(F(i,1),3)+V(F(i,2),3)+V(F(i,3),3))/3; % center of mass z-position
end
CoMface = [CoMxx' CoMyy' CoMzz']; % center of mass of each face

%% Step 3: Determine the centroid of the total surface
for i = 1:size(F,1)
    CoMface(i,:) = Aface(i).*CoMface(i,:);
end
Centroid = 1/sum(Aface).*sum(CoMface,1);