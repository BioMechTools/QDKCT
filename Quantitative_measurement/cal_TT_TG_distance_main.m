%% this file is to calculate TT-TG based on manual TG, and manual or auto reference line
clear all; clc; close all hidden;
% addpath('D:\project\4DTTTG\matlabCode\v0.5\coordinate\Source\');
%% functions location
addpath('.\function\')
[str_sourceNew, str_figure,str_figure_landmark,str_manual] = set_path();
%% Data Input

for numSubject = 1:8
    if numSubject>1
        n_subj = numSubject+1;
    else
        n_subj = numSubject;
    end
    str_subject = ['S00' num2str(n_subj)];
    for numSide = 1:2
        if numSide==1
            str_side = 'R';
        else
%             if n_subj == 4
%                 continue;
%             end
            str_side = 'L';
        end
        
        %%% read the manual TT and TG in static
        str_landmarks_file =  [str_manual 'Leo1\Subject' num2str(n_subj) '_' str_side '.txt' ] ;
        matLandmarks = load(str_landmarks_file);
        staticLandmarks1 = matLandmarks(1:4,1:3);
        str_landmarks_file =  [str_manual 'Leo2\Subject' num2str(n_subj) '_' str_side '.txt' ] ;
        matLandmarks = load(str_landmarks_file);
        staticLandmarks2 = matLandmarks(1:4,1:3);
        str_landmarks_file =  [str_manual 'Leo3\Subject' num2str(n_subj) '_' str_side '.txt' ] ;
        matLandmarks = load(str_landmarks_file);
        staticLandmarks3 = matLandmarks(1:4,1:3);
        str_landmarks_file =  [str_manual 'Seb\Subject' num2str(n_subj) '_' str_side '.txt' ] ;
        matLandmarks = load(str_landmarks_file);
        staticLandmarks4 = matLandmarks(1:4,1:3);
        str_landmarks_file =  [str_manual 'Mar\Subject' num2str(n_subj) '_' str_side '.txt' ] ;
        matLandmarks = load(str_landmarks_file);
        staticLandmarks5 = matLandmarks(1:4,1:3);       
        %%% stl
%         femPath = [str_source str_subject '\stlNew\Femur_' str_side '.stl' ] ;
%         patPath = [str_source str_subject '\stlNew\Patella_' str_side '.stl' ] ;
%         tibPath = [str_source str_subject '\stlNew\Tibia_' str_side '.stl' ] ;
%         [Tibia.F,Tibia.V,~] = fImportSTL(tibPath);
%         [Femur.F,Femur.V,~] = fImportSTL(femPath);
%         [Patella.F,Patella.V,~] = fImportSTL(patPath);
        %%% input coordinates
        fem_coord_path = [str_sourceNew str_subject '\matlab\Coordinate\StaticCoordinateFemur' str_side] ;
        tib_coord_path = [str_sourceNew str_subject '\matlab\Coordinate\StaticCoordinateTibia' str_side] ;
        load(fem_coord_path);
        load(tib_coord_path);
        femCoords.PD = StaticCoordinateFemur(2,1:3);
        femCoords.AP = StaticCoordinateFemur(1,1:3);
        femCoords.ML = StaticCoordinateFemur(3,1:3);
        femCoords.Origin = StaticCoordinateFemur(4,1:3);
        tibCoords.PD = StaticCoordinateTibia(2,1:3);
        tibCoords.AP = StaticCoordinateTibia(1,1:3);
        tibCoords.ML = StaticCoordinateTibia(3,1:3);
        tibCoords.Origin = StaticCoordinateTibia(4,1:3);
        
        % rotationMatrixTib = [tibCoords.PD tibCoords.ML tibCoords.AP];
        % getTT_TibiaNew(  R_axes,node,faceT )
        %% Reference line and TG on posterior condyle
        rotationMatrix = [femCoords.PD' femCoords.ML' femCoords.AP'];
        CutMatrix = [tibCoords.PD' tibCoords.ML' tibCoords.AP'];%[femCoords.PD femCoords.ML femCoords.AP];%
        rotationMatrixTib = [tibCoords.PD' tibCoords.ML' tibCoords.AP'];
        CT_frame.PD = [0 0 1];
        CT_frame.ML = [1 0 0];
        CT_frame.AP = [0 1 0];
        
%         [0 1 0;0 0 1;1 0 0];
        load([str_sourceNew str_subject '\matlab\Landmarks\28Tibia_Landmark_RefLine_' str_side]);
        load([str_sourceNew str_subject '\matlab\Landmarks\28Tibia_Landmark_TGPoint_' str_side]);
        load([str_sourceNew str_subject '\matlab\Landmarks\28Tibia_Landmark_TTPoint_' str_side]);
        
        TT_TG_RefPoints(1:2,:) = [Landmark_RefLine.LinePoint_Medial_3d';Landmark_RefLine.LinePoint_Lateral_3d'];
        TT_TG_RefPoints(3,:) = Landmark_TGPoint.TGPoint3D;
        TT_TG_RefPoints(4,:) = Landmark_TTPoint;
        
%         [TTTGdisList_Auto_Auto] = cal_TT_TG_dis(TT_TG_RefPoints(1:2,:),TT_TG_RefPoints(3,:),TT_TG_RefPoints(4,:),femCoords.ML,str_side)
        [TTTGdisList_Auto_Auto] = cal_TT_TG_dis(TT_TG_RefPoints(1:2,:),TT_TG_RefPoints(3,:),TT_TG_RefPoints(4,:),femCoords.ML,str_side)
        [Epsilon,theta,D,d] = get_para_Yao(TT_TG_RefPoints,tibCoords,tibCoords);
        [TT_TG_correctedAuto] = get_TT_TG_Yao(Epsilon,theta,D,d,TTTGdisList_Auto_Auto,str_side);
        %%% User1_1
        [TTTGdisList_Manual1_Manual1] = cal_TT_TG_dis(staticLandmarks1(1:2,:),staticLandmarks1(3,:),staticLandmarks1(4,:),femCoords.ML,str_side)
        [Epsilon,theta,D,d] = get_para_Yao(staticLandmarks1,CT_frame,tibCoords);
        [TT_TG_corrected1] = get_TT_TG_Yao(Epsilon,theta,D,d,TTTGdisList_Manual1_Manual1,str_side);
        %%% User1_2
        [TTTGdisList_Manual2_Manual2] = cal_TT_TG_dis(staticLandmarks2(1:2,:),staticLandmarks2(3,:),staticLandmarks2(4,:),femCoords.ML,str_side)
        [Epsilon,theta,D,d] = get_para_Yao(staticLandmarks2,CT_frame,tibCoords);
        [TT_TG_corrected2] = get_TT_TG_Yao(Epsilon,theta,D,d,TTTGdisList_Manual2_Manual2,str_side);
        %%% User1_3
        [TTTGdisList_Manual3_Manual3] = cal_TT_TG_dis(staticLandmarks3(1:2,:),staticLandmarks3(3,:),staticLandmarks3(4,:),femCoords.ML,str_side)
        [Epsilon,theta,D,d] = get_para_Yao(staticLandmarks3,CT_frame,tibCoords);
        [TT_TG_corrected3] = get_TT_TG_Yao(Epsilon,theta,D,d,TTTGdisList_Manual3_Manual3,str_side);
        %%% User2
        [TTTGdisList_Manual4_Manual4] = cal_TT_TG_dis(staticLandmarks4(1:2,:),staticLandmarks4(3,:),staticLandmarks4(4,:),femCoords.ML,str_side)
        [Epsilon,theta,D,d] = get_para_Yao(staticLandmarks4,CT_frame,tibCoords);
        [TT_TG_corrected4] = get_TT_TG_Yao(Epsilon,theta,D,d,TTTGdisList_Manual4_Manual4,str_side);
        %%% User3
        [TTTGdisList_Manual5_Manual5] = cal_TT_TG_dis(staticLandmarks5(1:2,:),staticLandmarks5(3,:),staticLandmarks5(4,:),femCoords.ML,str_side)
        [Epsilon,theta,D,d] = get_para_Yao(staticLandmarks5,CT_frame,tibCoords);
        [TT_TG_corrected5] = get_TT_TG_Yao(Epsilon,theta,D,d,TTTGdisList_Manual5_Manual5,str_side);  
        
        % str_file = [str_sourceNew 'TT_TG\' num2str(n_subj) 'TT_TG_dis' str_side];
        TT_TG_dis_list = [TTTGdisList_Auto_Auto,TTTGdisList_Manual1_Manual1,...
            TTTGdisList_Manual2_Manual2,TTTGdisList_Manual3_Manual3,...
            TTTGdisList_Manual4_Manual4,TTTGdisList_Manual5_Manual5,...
            TT_TG_correctedAuto,TT_TG_corrected1,TT_TG_corrected2,TT_TG_corrected3,...
            TT_TG_corrected4,TT_TG_corrected5,...
            Epsilon,theta,D,d];
%         TT_TG_dis_corrected_list = [TTTGdisList_Auto_Auto,];
        % save(str_file,'TT_TG_dis_list');
        if strcmp(str_side,'L')
            n_LR = 0;
        else
            n_LR = 1;
        end
        str_Column = char(66+(n_subj-1)*2+n_LR);
        filename = [str_manual 'TT_TG\TT_TG_tibia28.xlsx'];
        xlswrite(filename,TT_TG_dis_list',1,[str_Column '3:' str_Column '18']);
    end
end
