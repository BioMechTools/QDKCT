clear all; close all hidden; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This main function is to visualize the kinematics of DKCT data within
%   the points of tibial tubercle(TT), trochlear groove(TG), end points of
%   reference line(RL), the reference line and the projected line bewtween
%   TT-TG 
addpath('..\')
%% folder definition
str_sourceNew = '..\\..\\Data\\ProcessedData\\';
str_figure_landmark = '..\\..\\Data\\Draw\\landmarks\\';
str_figure = '..\\..\\Data\\Draw\\raw\\';
str_manual = '..\\..\\Data\\ManualTT_TG\\';

for n_subj = 7:7
    str_Subject = ['S00' num2str(n_subj)];
    for numSide = 1:2
        if numSide==1
            str_BoneSide = 'R';
        else
            str_BoneSide = 'L';
        end
        
        str_FolderCoordiate = [str_sourceNew str_Subject  '\matlab\Coordinate\'];% rotation matrix
        str_Folder_rotate = [str_sourceNew str_Subject '\matlab\Transform\'];

        %% Input of Coordiantes 
        load([str_FolderCoordiate 'StaticCoordinateFemur' str_BoneSide]);%Dynamic1FemurRT%
        load([str_FolderCoordiate 'StaticCoordinateTibia' str_BoneSide]);%Dynamic1PatellaRT%
        load([str_FolderCoordiate 'StaticCoordinatePatella' str_BoneSide]);%Dynamic1TibiaRT%

        load([str_Folder_rotate 'rigidD0toD1Femur' str_BoneSide]); 
        [regstr_matrixFem] = convertCPD2Affine3D(RigidS);
        Dynamic1CoordinateFemur = StaticCoordinateFemur;
        Dynamic1CoordinateFemur(4,1:3) = transformPointsForward(regstr_matrixFem,Dynamic1CoordinateFemur(4,1:3));
        Dynamic1CoordinateFemur(1:3,4) = Dynamic1CoordinateFemur(4,1:3)';
        Dynamic1CoordinateFemur(4,1:3) = 0;
        rotateFem = regstr_matrixFem;
        rotateFem.T(4,1:3) = 0; 
        Dynamic1CoordinateFemur(1:3,1:3) = transformPointsForward(rotateFem,Dynamic1CoordinateFemur(1:3,1:3));
        % temp = StaticCoordinatePatella(1:3,1:3);
        % StaticCoordinatePatella(1:3,1:3) = temp';

        load([str_Folder_rotate 'rigidD0toD1Patella' str_BoneSide]);  
        [regstr_matrixPat] = convertCPD2Affine3D(RigidS);
        Dynamic1CoordinatePatella = StaticCoordinatePatella;
        Dynamic1CoordinatePatella(4,1:3) = transformPointsForward(regstr_matrixPat,Dynamic1CoordinatePatella(4,1:3));
        Dynamic1CoordinatePatella(1:3,4) = Dynamic1CoordinatePatella(4,1:3)';
        Dynamic1CoordinatePatella(4,1:3) = 0;
        rotatePat = regstr_matrixPat;
        rotatePat.T(4,1:3) = 0; 
        Dynamic1CoordinatePatella(1:3,1:3) = transformPointsForward(rotatePat,Dynamic1CoordinatePatella(1:3,1:3));
        % temp = StaticCoordinateTibia(1:3,1:3);
        % StaticCoordinateTibia(1:3,1:3) = temp';

        load([str_Folder_rotate 'rigidD0toD1Tibia' str_BoneSide]);  
        [regstr_matrixTib] = convertCPD2Affine3D(RigidS);
        Dynamic1CoordinateTibia = StaticCoordinateTibia;
        Dynamic1CoordinateTibia(4,1:3) = transformPointsForward(regstr_matrixTib,Dynamic1CoordinateTibia(4,1:3));
        Dynamic1CoordinateTibia(1:3,4) = Dynamic1CoordinateTibia(4,1:3)';
        Dynamic1CoordinateTibia(4,1:3) = 0;
        rotateTib = regstr_matrixTib;
        rotateTib.T(4,1:3) = 0; 
        Dynamic1CoordinateTibia(1:3,1:3) = transformPointsForward(rotateTib,Dynamic1CoordinateTibia(1:3,1:3));
        save([str_FolderCoordiate 'Dynamic1CoordinateFemur' str_BoneSide], 'Dynamic1CoordinateFemur');
        save([str_FolderCoordiate 'Dynamic1CoordinatePatella' str_BoneSide], 'Dynamic1CoordinatePatella');
        save([str_FolderCoordiate 'Dynamic1CoordinateTibia' str_BoneSide], 'Dynamic1CoordinateTibia');
    end
end