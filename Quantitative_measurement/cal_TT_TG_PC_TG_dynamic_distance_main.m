%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This main function is to calculate the kinematics of DKCT data

clear all; clc;close all;
%% files location
addpath('..\visualization\')
% addpath('..\..\')
% addpath('..\..\Kinematics_Cal\')
% addpath('..\..\function\')
% str_source = 'I:\UT\DKCT\ProcessDataNew\';
addpath('..\')
%% folder definition
str_sourceNew = '..\\..\\Data\\ProcessedData\\';
str_figure_landmark = '..\\..\\Data\\Draw\\landmarks\\';
str_figure = '..\\..\\Data\\Draw\\raw\\';
str_manual = '..\\..\\Data\\ManualTT_TG\\';

for n_subj = 7:7
    str_subject = ['S00' num2str(n_subj)];
    for numSide = 1:2
        if numSide==1
            str_BoneSide = 'R';
        else
            str_BoneSide = 'L';
        end
        str_Dynamic2DynamicFolder = [str_sourceNew str_subject '\\matlab\\Transform\\'];
        str_FolderCoordiate = [str_sourceNew str_subject  '\matlab\Coordinate\'];
        str_TT_TG_Folder = [str_sourceNew str_subject '\matlab\TT_TG\'];% rotation matrix
        %% Input of stl        
        %%% stl
        femPath = [str_sourceNew str_subject '\stlNew\DynamicResample\Femur_' str_BoneSide '.stl' ] ;
        patPath = [str_sourceNew str_subject '\stlNew\DynamicResample\Patella_' str_BoneSide '.stl' ] ;
        tibPath = [str_sourceNew str_subject '\stlNew\DynamicResample\Tibia_' str_BoneSide '.stl' ] ;
        [Tibia.F,Tibia.V,~] = fImportSTL(tibPath);
        [Femur.F,Femur.V,~] = fImportSTL(femPath);
        [Patella.F,Patella.V,~] = fImportSTL(patPath);
        %% Input of Coordiantes
        %%% the structure of the femur coordinate is 
        %%% [referenceFrames.f.X;
        %%%  referenceFrames.f.Y;
        %%%  referenceFrames.f.Z;
        %%%  referenceFrames.f.origin]
        load([str_FolderCoordiate 'StaticCoordinateFemur' str_BoneSide]);%Dynamic1FemurRT%
        load([str_FolderCoordiate 'StaticCoordinateTibia' str_BoneSide]);%Dynamic1PatellaRT%
        load([str_FolderCoordiate 'StaticCoordinatePatella' str_BoneSide]);%Dynamic1TibiaRT%
        %% Femur coordinate (this is aimed to adjust the sign of TT-TG distance)
%         [femCoords]=cal_Femur_direction(Femur, Tibia, Patella);
        femCoords.PD = StaticCoordinateFemur(2,1:3);
        femCoords.AP = StaticCoordinateFemur(1,1:3);
        femCoords.ML = StaticCoordinateFemur(3,1:3);
        femCoords.Origin = StaticCoordinateFemur(4,1:3);
        tibCoords.PD = StaticCoordinateTibia(2,1:3);
        tibCoords.AP = StaticCoordinateTibia(1,1:3);
        tibCoords.ML = StaticCoordinateTibia(3,1:3);
        tibCoords.Origin = StaticCoordinateTibia(4,1:3);
        %% Kinematics calculation
        %%% Femur & Tibia
        A= repmat(eye(4,4), 1, 1, 12); B=A;
        A(1:3,1:4,1)= [StaticCoordinateFemur(1:3, 1:3)' StaticCoordinateFemur(4,1:3)'];
        B(1:3,1:4,1)= [StaticCoordinateTibia(1:3, 1:3)' StaticCoordinateTibia(4,1:3)'];
        currentCoordinateFemur = StaticCoordinateFemur;
        currentCoordinateTibia = StaticCoordinateTibia;
        currentCoordinatePatella = StaticCoordinatePatella;
        %% read the landmarks
        %%% read the landmarks
        str_TG_file = [str_sourceNew str_subject '\matlab\Landmarks\28Tibia_Landmark_TGPoint_' str_BoneSide];
        str_PCL_file = [str_sourceNew str_subject '\matlab\Landmarks\28Tibia_Landmark_RefLine_' str_BoneSide];
        str_TT_file = [str_sourceNew str_subject '\matlab\Landmarks\28Tibia_Landmark_TTPoint_' str_BoneSide];
        load(str_TG_file);
        load(str_PCL_file);
        load(str_TT_file);
        TT_TG_RefPoints = [Landmark_RefLine.LinePoint_Medial_3d';Landmark_RefLine.LinePoint_Lateral_3d';Landmark_TGPoint.TGPoint3D';Landmark_TTPoint];
        TT_TG_RefPointsOrigin = TT_TG_RefPoints;
        TTTGdisList = zeros(12,1);
        TTTGdisList(1)=cal_TT_TG_dis(TT_TG_RefPoints(1:2,:),TT_TG_RefPoints(3,:),TT_TG_RefPoints(4,:),femCoords.ML,str_BoneSide);
        PatTGdisList = zeros(12,1);
        PatTGdisList(1)=cal_TT_TG_dis(TT_TG_RefPoints(1:2,:),TT_TG_RefPoints(3,:),currentCoordinatePatella(4,1:3),femCoords.ML,str_BoneSide);
        drawSTL_TT_TG( Femur, Tibia, Patella, ...
            currentCoordinatePatella(4,1:3),TT_TG_RefPoints(3,:),TT_TG_RefPoints(1:2,:))

        for num =0:10
            %%% load rotation matrix and convert to affine3d
            str_Frame2Frame_femur = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dFemur' str_BoneSide],num,num+1);
            load(str_Frame2Frame_femur);
            [regstr_matrixFem] = convertCPD2Affine3D(RigidS);
            str_Frame2Frame_tibia = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dTibia' str_BoneSide],num,num+1);
            load(str_Frame2Frame_tibia);
            [regstr_matrixTib] = convertCPD2Affine3D(RigidS);
            str_Frame2Frame_patella = sprintf([str_Dynamic2DynamicFolder 'rigidD%dtoD%dPatella' str_BoneSide],num,num+1);
            load(str_Frame2Frame_patella);
            [regstr_matrixPat] = convertCPD2Affine3D(RigidS);

            [currentCoordinateFemur] = transform_dynamic2dynamic(currentCoordinateFemur,regstr_matrixFem);
            [currentCoordinateTibia] = transform_dynamic2dynamic(currentCoordinateTibia,regstr_matrixTib); 
        %     RotateFem = regstr_matrixFem.T;
        %     RotateTib = regstr_matrixTib.T;
            A(1:3,1:4,num + 2)= [currentCoordinateFemur(1:3,1:3)' currentCoordinateFemur(4,1:3)'];
            B(1:3,1:4,num + 2)= [currentCoordinateTibia(1:3,1:3)' currentCoordinateTibia(4,1:3)'];

            %%%TTTG
            TT_TG_RefPoints(1:3,:)=transformPointsForward(regstr_matrixFem,TT_TG_RefPoints(1:3,:));
            TT_TG_RefPoints(4,:)=transformPointsForward(regstr_matrixTib,TT_TG_RefPoints(4,:));
            TTTGdisList(num+2)=cal_TT_TG_dis(TT_TG_RefPoints(1:2,:),TT_TG_RefPoints(3,:),TT_TG_RefPoints(4,:),femCoords.ML,str_BoneSide);
            %%%PatTG
            currentCoordinatePatella(4,1:3)=transformPointsForward(regstr_matrixPat,currentCoordinatePatella(4,1:3));
            PatTGdisList(num+2)=cal_TT_TG_dis(TT_TG_RefPoints(1:2,:),TT_TG_RefPoints(3,:),currentCoordinatePatella(4,1:3),femCoords.ML,str_BoneSide);

            drawSTL_TT_TG( Femur, Tibia, Patella, ...
            currentCoordinatePatella(4,1:3),TT_TG_RefPoints(3,:),TT_TG_RefPoints(1:2,:));
        end
        [TFrot, TFtrans]= calcJCSkin(A, B, lower(str_BoneSide));
%         if strcmp(str_BoneSide,'L')
%             PatTGdisList = -PatTGdisList;
%             TTTGdisList = -TTTGdisList;
%         end
        save([str_TT_TG_Folder 'TT_TG_distance_List_' str_BoneSide],'TTTGdisList');
        save([str_TT_TG_Folder 'PC_TG_distance_List_' str_BoneSide],'PatTGdisList');
        save([str_TT_TG_Folder 'TFrot_List_' str_BoneSide],'TFrot');
    end
end
