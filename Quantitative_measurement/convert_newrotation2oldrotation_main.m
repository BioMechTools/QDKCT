%% this file is aimed to convert the new roatation matrix (S to 5) to the original rotation matrix (S to 1)
clear;clc;close all hidden;
addpath('..\')
addpath('..\basic_function\')
[str_input_folder, str_figure,str_figure_landmark,str_manual] = set_path();
%% folder
n_subj = 7;
% str_source = 'I:\UT\DKCT\ProcessData\';
str_BoneSide = 'L';
str_Subject = ['S00' num2str(n_subj)];
% str_input_folder = 'I:\\UT\\DKCT\\ProcessDataNew\\';
str_rotation_folder = [str_input_folder str_Subject '\\matlab\\TransformPreNew\\'];

str_output_folder = [str_input_folder str_Subject '\\matlab\\Transform\\'];
%% convert the static to D5 to static to D1
str_static2D5 = [str_rotation_folder 'rigidStatictoD5Femur' str_BoneSide];
load(str_static2D5);
[rot_Femur] = convertCPD2Affine3D(RigidS);
str_static2D5 = [str_rotation_folder 'rigidStatictoD5Tibia' str_BoneSide];
load(str_static2D5);
[rot_Tibia] = convertCPD2Affine3D(RigidS);
str_static2D5 = [str_rotation_folder 'rigidStatictoD5Patella' str_BoneSide];
load(str_static2D5);
[rot_Patella] = convertCPD2Affine3D(RigidS);

for num = 5:-1:2
    str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dFemur' str_BoneSide],num,num-1);
    load(str_Frame2Frame);
    [rot_temp] = convertCPD2Affine3D(RigidS);
    rot_Femur.T = rot_Femur.T*rot_temp.T;
    str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dTibia' str_BoneSide],num,num-1);
    load(str_Frame2Frame);
    [rot_temp] = convertCPD2Affine3D(RigidS);
    rot_Tibia.T = rot_Tibia.T*rot_temp.T;
    str_Frame2Frame = sprintf([str_rotation_folder 'rigidD%dtoD%dPatella' str_BoneSide],num,num-1);
    load(str_Frame2Frame);
    [rot_temp] = convertCPD2Affine3D(RigidS);
    rot_Patella.T = rot_Patella.T*rot_temp.T;
end
%% test
%%% stl
% femPath = [str_source str_Subject '\stlNew\Femur_' str_BoneSide '.stl' ] ;
% patPath = [str_source str_Subject '\stlNew\Patella_' str_BoneSide '.stl' ] ;
% tibPath = [str_source str_Subject '\stlNew\Tibia_' str_BoneSide '.stl' ] ;
% [Tibia.faces,Tibia.vertices,~] = fImportSTL(tibPath);
% [Femur.faces,Femur.vertices,~] = fImportSTL(femPath);
% [Patella.faces,Patella.vertices,~] = fImportSTL(patPath);
% [Femur.vertices] = transformPointsForward(rot_Femur,Femur.vertices); 
% [Tibia.vertices] = transformPointsForward(rot_Tibia,Tibia.vertices);
% [Patella.vertices] = transformPointsForward(rot_Patella,Patella.vertices);
% 
% stlwrite('.\femur.stl',Femur);
% stlwrite('.\tibia.stl',Tibia);
% stlwrite('.\patella.stl',Patella);
%% save the rotation matrix from static to D1
[RigidS] = convertAffine3D2CPD(rot_Femur);
str_output = sprintf([str_output_folder 'rigidD0toD1Femur' str_BoneSide]);
save(str_output,'RigidS');

[RigidS] = convertAffine3D2CPD(rot_Tibia);
str_output = sprintf([str_output_folder 'rigidD0toD1Tibia' str_BoneSide]);
save(str_output,'RigidS');

[RigidS] = convertAffine3D2CPD(rot_Patella);
str_output = sprintf([str_output_folder 'rigidD0toD1Patella' str_BoneSide]);
save(str_output,'RigidS');


%% convert the D2toD1 to D1toD2
for num = 1:4
    str_Frame2Frame_femur = sprintf([str_rotation_folder 'rigidD%dtoD%dFemur' str_BoneSide],num+1,num);
    load(str_Frame2Frame_femur);
    [rot_temp] = convertCPD2Affine3D(RigidS);
    temp = inv(rot_temp.T);
    temp(1:3,4)=0;
    temp(4,4)=1;
    rot_temp.T = temp;
    [RigidS] = convertAffine3D2CPD(rot_temp);
    str_output = sprintf([str_output_folder 'rigidD%dtoD%dFemur' str_BoneSide],num,num+1);
    save(str_output,'RigidS');
    
    str_Frame2Frame_femur = sprintf([str_rotation_folder 'rigidD%dtoD%dTibia' str_BoneSide],num+1,num);
    load(str_Frame2Frame_femur);
    [rot_temp] = convertCPD2Affine3D(RigidS);
    temp = inv(rot_temp.T);
    temp(1:3,4)=0;
    temp(4,4)=1;
    rot_temp.T = temp;
    [RigidS] = convertAffine3D2CPD(rot_temp);
    str_output = sprintf([str_output_folder 'rigidD%dtoD%dTibia' str_BoneSide],num,num+1);
    save(str_output,'RigidS');
    
    str_Frame2Frame_femur = sprintf([str_rotation_folder 'rigidD%dtoD%dPatella' str_BoneSide],num+1,num);
    load(str_Frame2Frame_femur);
    [rot_temp] = convertCPD2Affine3D(RigidS);
    temp = inv(rot_temp.T);
    temp(1:3,4)=0;
    temp(4,4)=1;
    rot_temp.T = temp;
    [RigidS] = convertAffine3D2CPD(rot_temp);
    str_output = sprintf([str_output_folder 'rigidD%dtoD%dPatella' str_BoneSide],num,num+1);
    save(str_output,'RigidS');
end
