function [refL,VcsL,IA] = computeReferenceLength(V)
%Computes the reference length from the vertices of a contour around the
%condyles. The reference length is the mean side length of the minimal
%bounding box. The condyles are found by finding the location of maximal
%crossecitonal area along the first principle component.
%
%   [refL] = COMPUTEREFERENCELENGTH(ForVwide,V,gridSize)
%
%--Input
%   V           -    Vertices of the distal femur
%
%--Output
%   refL        -   Reference length
%   VcsL        -   Vertices in pca space
%   IA          -   pca's 'inertial' axes
%
%  (V1.1) ORL Nijmegen, Max Bakker 2016

%Changelog:
%       performance increase and more sensible IO

%now, with pca in stead of inertial axes. saves ~2s
IA = pca(V) ;
%[~,IA,~] = computeInertialAxes(ForVwide,V,gridSize) ;
[CA,sectionCenters,~,VcsL] = computeCrossSectionAreaAlongAxis(V,IA,2) ;
[~,peakI] = max(CA) ;
VwideI = abs(VcsL(:,1) - sectionCenters(peakI)*range(VcsL(:,1))/100 - min(VcsL(:,1) )) < 1 ;
Vcontour = VcsL(VwideI,[2,3]);


boxSideL = findBoundingBox2D(Vcontour,@prod,@min )  ;
refL = mean(boxSideL) ;
end