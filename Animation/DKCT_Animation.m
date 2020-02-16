%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This main function is to visualize the kinematics of DKCT data within
%   the points of tibial tubercle(TT), trochlear groove(TG), end points of
%   reference line(RL), the reference line and the projected line bewtween
%   TT-TG 
clear all; clc;close all;
%% files location
addpath('..\')
addpath('..\visualization\')
%% folder definition
str_sourceNew = '..\\..\\Data\\ProcessedData\\';
str_figure_landmark = '..\\..\\Data\\Draw\\landmarks\\';
str_figure = '..\\..\\Data\\Draw\\raw\\';
str_manual = '..\\..\\Data\\ManualTT_TG\\';

str_BoneSide = 'R';
n_subj = 7;
str_Subject = ['S00' num2str(n_subj)];
str_Dynamic2DynamicFolder = [str_sourceNew str_Subject '\\matlab\\Transform\\'];
%points position of TTTG/6x3 matrix (point 1 of RL; point 2 of RL;TG on Fem;TG on bottom of Pat,TG on middle of Pat; TT on Tib )
% str_sourceTTTG = 'D:\project\4DTTTG\data\Output\1\';
%% Input of point of TT,TG and end points of reference line
%%% the position of the points should be the position of static scan
load([str_sourceNew str_Subject '\matlab\Landmarks\28Tibia_Landmark_RefLine_' str_BoneSide]);
load([str_sourceNew str_Subject '\matlab\Landmarks\28Tibia_Landmark_TGPoint_' str_BoneSide]);
load([str_sourceNew str_Subject '\matlab\Landmarks\28Tibia_Landmark_TTPoint_' str_BoneSide]);

TT_TG_RefPoints(1:2,:) = [Landmark_RefLine.LinePoint_Medial_3d';Landmark_RefLine.LinePoint_Lateral_3d'];
TT_TG_RefPoints(3,:) = Landmark_TGPoint.TGPoint3D;
TT_TG_RefPoints(4,:) = Landmark_TTPoint;
        
TT = TT_TG_RefPoints(4,:);
TG_Femur = TT_TG_RefPoints(3,:);
pointRL = TT_TG_RefPoints(1:2,:);
TG_PatM = TT_TG_RefPoints(4,:);         
TG_PatB = TT_TG_RefPoints(4,:);   

%% Input of STL 
femPath = [str_sourceNew str_Subject '\stlNew\DynamicResample\Femur_' str_BoneSide '_Static_cut.stl' ] ;
patPath = [str_sourceNew str_Subject '\stlNew\DynamicResample\Patella_' str_BoneSide '_Static.stl' ] ;
tibPath = [str_sourceNew str_Subject '\stlNew\DynamicResample\Tibia_' str_BoneSide '_Static_cut.stl' ] ;
FemurSTL= stlread(femPath);
PatellaSTL= stlread(patPath);
TibiaSTL= stlread(tibPath);
%%% transform the stl in static position to the first frame of DKCT scan
str_Frame2Frame_femur = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dFemur' str_BoneSide],0,1);
load(str_Frame2Frame_femur);
[regstr_matrixFem] = convertCPD2Affine3D(RigidS);
% load([str_sourceR '00MatrixFemur']); 
% load([str_sourceR '00MatrixTibia']);  
% load([str_sourceR '00MatrixPatella']);  
str_Frame2Frame_tibia = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dTibia' str_BoneSide],0,1);
load(str_Frame2Frame_tibia);
[regstr_matrixTib] = convertCPD2Affine3D(RigidS);
str_Frame2Frame_patella = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dPatella' str_BoneSide],0,1);
load(str_Frame2Frame_patella);
[regstr_matrixPat] = convertCPD2Affine3D(RigidS);

FemurSTL.vertices = transformPointsForward(regstr_matrixFem,FemurSTL.vertices);
PatellaSTL.vertices = transformPointsForward(regstr_matrixPat,PatellaSTL.vertices);
TibiaSTL.vertices = transformPointsForward(regstr_matrixTib,TibiaSTL.vertices);

%% TT TG RL from static to dynamic
TG = transformPointsForward(regstr_matrixFem,TG_Femur);
pointRL(1,:) = transformPointsForward(regstr_matrixFem,pointRL(1,:));
pointRL(2,:) = transformPointsForward(regstr_matrixFem,pointRL(2,:));
TT = transformPointsForward(regstr_matrixTib,TT);
TG_PatM = transformPointsForward(regstr_matrixTib,TG_PatM);
TG_PatB = transformPointsForward(regstr_matrixTib,TG_PatB);

%% Input of Coordiantes 
load([str_sourceNew str_Subject '\matlab\Coordinate\Dynamic1CoordinateFemur' str_BoneSide]);%Dynamic1FemurRT%
load([str_sourceNew str_Subject '\matlab\Coordinate\Dynamic1CoordinateTibia' str_BoneSide]);%Dynamic1PatellaRT%
load([str_sourceNew str_Subject '\matlab\Coordinate\Dynamic1CoordinatePatella' str_BoneSide']);%Dynamic1TibiaRT%
temp = Dynamic1CoordinateFemur(1:3,1:3);
Dynamic1CoordinateFemur(1:3,1:3) = temp';
FemurCoordsPosDir = Dynamic1CoordinateFemur;

temp = Dynamic1CoordinatePatella(1:3,1:3);
Dynamic1CoordinatePatella(1:3,1:3) = temp';
PatellaCoordsPosDir = Dynamic1CoordinatePatella;

temp = Dynamic1CoordinateTibia(1:3,1:3);
Dynamic1CoordinateTibia(1:3,1:3) = temp';
TibiaCoordsPosDir = Dynamic1CoordinateTibia;

%% Input of rotation matrix
FemurRMs = zeros(4,4,11);
PatellaRMs = zeros(4,4,11);
TibiaRMs = zeros(4,4,11);
% the first transformation is itself
FemurRMs(:,:,1) = diag([1 1 1 1]);
PatellaRMs(:,:,1) = diag([1 1 1 1]);
TibiaRMs(:,:,1) = diag([1 1 1 1]);
% str_OutputTest = sprintf(['G:\\Hao\\DKCT\\ProcessData\\S001\\stl\\Test\\' 'Femur%d.stl'],1);
% stlwrite(str_OutputTest,FemurSTL);
% str_OutputTest = sprintf(['G:\\Hao\\DKCT\\ProcessData\\S001\\stl\\Test\\' 'Patella%d.stl'],1);
% stlwrite(str_OutputTest,PatellaSTL);
% str_OutputTest = sprintf(['G:\\Hao\\DKCT\\ProcessData\\S001\\stl\\Test\\' 'Tibia%d.stl'],1);
% stlwrite(str_OutputTest,TibiaSTL);
for num =1:10
    %%% femur
    str_Frame2Frame_femur = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dFemur' str_BoneSide],num,num+1);
    load(str_Frame2Frame_femur);
    [regstr_matrixFem] = convertCPD2Affine3D(RigidS);
    FemurRMs(:,:,num+1) = FemurRMs(:,:,num)*regstr_matrixFem.T;
    %%% patella
    str_Frame2Frame_patella = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dPatella' str_BoneSide],num,num+1);
    load(str_Frame2Frame_patella);
    [regstr_matrixPat] = convertCPD2Affine3D(RigidS);
%     load([str_sourceR sprintf('%02dMatrixPatella',num+2)]);
    PatellaRMs(:,:,num+1) = PatellaRMs(:,:,num)*(regstr_matrixPat.T);
    %%% tibia
    str_Frame2Frame_tibia = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dTibia' str_BoneSide],num,num+1);
    load(str_Frame2Frame_tibia);
    [regstr_matrixTib] = convertCPD2Affine3D(RigidS);
%     load([str_sourceR sprintf('%02dMatrixTibia',num+2)]);
    TibiaRMs(:,:,num+1) = TibiaRMs(:,:,num)*(regstr_matrixTib.T);
    
end

%% Draw

% TTNew = transformPointsInverse(regstr_matrixTib,TT_TG_RefPoints(6,:));
drawKinematics( FemurSTL, TibiaSTL, PatellaSTL, FemurRMs, TibiaRMs, PatellaRMs, ...
    FemurCoordsPosDir,TibiaCoordsPosDir,PatellaCoordsPosDir, TT, TG,TG_PatM,TG_PatB, pointRL,[],1);

