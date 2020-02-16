%% this file is aimed to transform obtain the transform matrix from dynamic frame to static

clear; close all hidden; clc;

addpath('..\')
addpath('..\basic_function\')
[str_input_folder, str_figure,str_figure_landmark,str_manual] = set_path();

for numSubject = 6:6
    if numSubject>1
        n_subj = numSubject+1;
    else
        n_subj = numSubject;
    end
    str_subject = ['S00' num2str(n_subj)];
    for numSide = 1:2
        if numSide==1
            str_BoneSide = 'R';
        else
            if n_subj == 4
                continue;
            end
            str_BoneSide = 'L';
        end
        str_Subject = ['S00' num2str(n_subj)];
        str_rotation_folder = [str_input_folder str_Subject '\\matlab\\TransformPreNew\\'];
        str_static2dynamic_folder = [str_input_folder str_Subject '\\matlab\\TransformStatic2Dynamic\\'];
        str_output_folder = [str_input_folder str_Subject '\\matlab\\Transform\\'];
        
        %%% input the static to dynamic 5
        str_static2D5 = [str_rotation_folder 'rigidStatictoD5Femur' str_BoneSide];
        load(str_static2D5);
        [rot_FemurO] = convertCPD2Affine3D(RigidS);
        str_static2D5 = [str_rotation_folder 'rigidStatictoD5Tibia' str_BoneSide];
        load(str_static2D5);
        [rot_TibiaO] = convertCPD2Affine3D(RigidS);
        str_static2D5 = [str_rotation_folder 'rigidStatictoD5Patella' str_BoneSide];
        load(str_static2D5);
        [rot_PatellaO] = convertCPD2Affine3D(RigidS);
        str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Femur_' str_BoneSide ],5);
        rot_Femur = rot_FemurO;
        rot_Tibia = rot_TibiaO;
        rot_Patella = rot_PatellaO;
        rot_S2D = rot_Femur;
        save(str_static2dynamic,'rot_S2D');
        str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Tibia_' str_BoneSide ],5);
        rot_S2D = rot_Tibia;
        save(str_static2dynamic,'rot_S2D');
        str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Patella_' str_BoneSide ],5);
        rot_S2D = rot_Patella;
        save(str_static2dynamic,'rot_S2D');
        
        for num = 5:-1:2
            str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dFemur' str_BoneSide],num,num-1);
            load(str_Frame2Frame);
            [rot_temp] = convertCPD2Affine3D(RigidS);
            rot_Femur.T = rot_Femur.T*rot_temp.T;
            str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Femur_' str_BoneSide ],num-1);
            rot_S2D = rot_Femur;
            save(str_static2dynamic,'rot_S2D');
            
            str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dTibia' str_BoneSide],num,num-1);
            load(str_Frame2Frame);
            [rot_temp] = convertCPD2Affine3D(RigidS);
            rot_Tibia.T = rot_Tibia.T*rot_temp.T;
            str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Tibia_' str_BoneSide ],num-1);
            rot_S2D = rot_Tibia;
            save(str_static2dynamic,'rot_S2D');
            
            str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dPatella' str_BoneSide],num,num-1);
            load(str_Frame2Frame);
            [rot_temp] = convertCPD2Affine3D(RigidS);
            rot_Patella.T = rot_Patella.T*rot_temp.T;
            str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Patella_' str_BoneSide ],num-1);
            rot_S2D = rot_Patella;
            save(str_static2dynamic,'rot_S2D');
        end
        rot_Femur = rot_FemurO;
        rot_Tibia = rot_TibiaO;
        rot_Patella = rot_PatellaO;     
        for num = 5:10
            str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dFemur' str_BoneSide],num,num+1);
            load(str_Frame2Frame);
            [rot_temp] = convertCPD2Affine3D(RigidS);
            rot_Femur.T = rot_Femur.T*rot_temp.T;
            str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Femur_' str_BoneSide ],num+1);
            rot_S2D = rot_Femur;
            save(str_static2dynamic,'rot_S2D');
            
            str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dTibia' str_BoneSide],num,num+1);
            load(str_Frame2Frame);
            [rot_temp] = convertCPD2Affine3D(RigidS);
            rot_Tibia.T = rot_Tibia.T*rot_temp.T;
            str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Tibia_' str_BoneSide ],num+1);
            rot_S2D = rot_Tibia;
            save(str_static2dynamic,'rot_S2D');
            
            str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dPatella' str_BoneSide],num,num+1);
            load(str_Frame2Frame);
            [rot_temp] = convertCPD2Affine3D(RigidS);
            rot_Patella.T = rot_Patella.T*rot_temp.T;
            str_static2dynamic = sprintf([str_static2dynamic_folder 'Static2D%d_Patella_' str_BoneSide ],num+1);
            rot_S2D = rot_Patella;
            save(str_static2dynamic,'rot_S2D');
        end
    end
end