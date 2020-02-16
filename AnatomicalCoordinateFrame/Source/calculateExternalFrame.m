function [ XYZest ] = calculateExternalFrame(femV,patV,tibV) 
%CALCULATEEXTERNALFRAME Estimate AP(X), PD(Y) and ML(Z) and the basis of 
%the relative locations of the femur, patella and tibia.
%   The underlying logic is that
%       The tibia is more distal than the femur
%       The patella is more anterior than both
%
%--Input
%   femV   -   Vertices of the femur
%   tibV   -   Vertices of the tibia
%   patV   -   Vertices of the patella
%       Note that it is essential that the three vertex sets are registred
%       to the same global space.
%
%--Output
%   XYZest -   struct with three field, having the coordinates of the
%   respective axis
%
%   see also ERCrefFrame

CoMfem = mean(femV) ;
CoMtib = mean(tibV) ;
CoMpat = mean(patV) ;

femTibVec = CoMtib - CoMfem;
femTibVec = femTibVec / rssq(femTibVec) ; %now NORMALIZE!
femPatVec = CoMpat-CoMfem ;
femPatVec = femPatVec / rssq(femPatVec) ;

XYZest.Y = -femTibVec ;

%This is the only slightly fancy thing here:
    %we determine AP(X) as the line connecting the patellar CM to the line
    %connecting the tibia and the femur(femTibVec) in a right angle 
XYZest.X = femPatVec + dot(femPatVec,femTibVec,1) ;
XYZest.X = XYZest.X / rssq(XYZest.X) ;
XYZest.Z = cross(XYZest.X,XYZest.Y) ;
XYZest.Z = XYZest.Z  / rssq(XYZest.Z ) ;

XYZest.X = cross(XYZest.Y,XYZest.Z) ;

end

