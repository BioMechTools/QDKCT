%This file is used to calculate the anatomical coordinates based on STL
%files
clear;close all hidden;clc
%% get the locations of the files
%add the functions to matlabs path 
addpath(['..' filesep 'coordinate' filesep 'Source']) ;
%also ad geom3d, used for some measures, but mostly for visualisation
addpath(['..' filesep 'coordinate' filesep 'geom3d']) ;
%filesep is either /  or \, depending on the OS
%% set the data path
str_main_path = 'I:\UT\DKCT_Quantitative_Measurement\Data\';
%% Target folder for stl
str_BoneSide = 'R';
str_Subject = 'S007';

femPath = [str_main_path str_Subject '\stlNew\DynamicResample\Femur_' str_BoneSide '_Static.stl' ] ;
patPath = [str_main_path str_Subject '\stlNew\DynamicResample\Patella_' str_BoneSide '_Static.stl' ] ;
tibPath = [str_main_path str_Subject '\stlNew\DynamicResample\Tibia_' str_BoneSide '_Static.stl'  ] ;

%% Do the calculation
[referenceFrames,stlData,DiagInfos] = ERCkneeReferenceFrames(femPath,patPath,tibPath) ;

%% save 
StaticCoordinateFemur = eye(4,4);
StaticCoordinateFemur(1:4,1:3) = [referenceFrames.f.X;referenceFrames.f.Y;referenceFrames.f.Z;referenceFrames.f.origin];
save([str_main_path str_Subject '\matlab\Coordinate\StaticCoordinateFemur' str_BoneSide '.mat'] ,'StaticCoordinateFemur');
StaticCoordinateTibia = eye(4,4);
StaticCoordinateTibia(1:4,1:3) = [referenceFrames.t.X;referenceFrames.t.Y;referenceFrames.t.Z;referenceFrames.t.origin];
save([str_main_path str_Subject '\matlab\Coordinate\StaticCoordinateTibia' str_BoneSide '.mat'] ,'StaticCoordinateTibia');
StaticCoordinatePatella = eye(4,4);
StaticCoordinatePatella(1:4,1:3) = [referenceFrames.p.X;referenceFrames.p.Y;referenceFrames.p.Z;referenceFrames.p.origin];
save([str_main_path str_Subject '\matlab\Coordinate\StaticCoordinatePatella' str_BoneSide '.mat'] ,'StaticCoordinatePatella');
save([str_main_path str_Subject '\matlab\Coordinate\DiagInfos'] ,'DiagInfos');

%% Below is the codes for visualization
%plot the results in the console
disp('Femoral reference Frame');
disp(referenceFrames.f)
disp('Patellar reference Frame');
disp(referenceFrames.p)
disp('Tibial reference Frame');
disp(referenceFrames.t)

%% Visualise the result
figure(1)
%plotPlaces handles the creation of axes in figures (it calls subplot)
hold on; axis equal ; axis tight ; view(3) ;
patch('faces',stlData(1).F,'Vertices',stlData(1).V,'edgeAlpha',0.5,'FaceAlpha',0.3);
plotCoords(referenceFrames.f)

patch('faces',stlData(2).F,'Vertices',stlData(2).V,'edgeAlpha',0.5,'FaceAlpha',0.3);
plotCoords(referenceFrames.p)

patch('faces',stlData(3).F,'Vertices',stlData(3).V,'edgeAlpha',0.5,'FaceAlpha',0.3);
plotCoords(referenceFrames.t)

%% visualise the steps. This is especially usefull with unexpected results

plotFemurDiagInfo(DiagInfos.f,stlData(1))
plotPatellaDiagInfo(DiagInfos.p,stlData(2)) ;
plotTibiaDiagInfo(DiagInfos.t,stlData(3)) ;

