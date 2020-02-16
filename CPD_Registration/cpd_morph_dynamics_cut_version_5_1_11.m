clear all; close all; clc;
addpath(['..' filesep 'CPD_registration' filesep 'function']) ;
addpath(['..' filesep 'CPD_registration' filesep 'core']) ;
str_folder = '..\\..\\Data\\';
str_DynamicFolder = '..\\..\\Data\\';
SubjectNumber = 'S007';
nFrameNum = 5;
str_BoneType = 'Patella';%%% choose patella, tibia, and femur
str_BoneSide = 'R';
outliers = 0.001;
fgt = 0;
for nFrameNum = 5:-1:2
    str_Target = sprintf('D%dtoD%d',nFrameNum,nFrameNum-1);
    str_subFolderRigid = [str_folder SubjectNumber '\\matlab\\TransformPreNew\\'];
    str_RigidFileOutput = sprintf([str_subFolderRigid 'rigid' str_Target str_BoneType str_BoneSide '.mat']);
    
    str_TargetSTL = sprintf([str_DynamicFolder SubjectNumber '\\stlNew\\DynamicResample\\' str_BoneType str_BoneSide '_'  num2str(nFrameNum-1) '.stl']);
    str_DistanceVOutput = sprintf([str_subFolderRigid 'Dist_V' str_Target str_BoneType str_BoneSide '.mat']);
    str_DistancesurfOutput = sprintf([str_subFolderRigid 'Dist_surf' str_Target str_BoneType str_BoneSide '.mat']);
    str_DistanceCorrespondenceOutput = sprintf([str_subFolderRigid 'Correspondence' str_Target str_BoneType str_BoneSide '.mat']);
    str_DistanceVerticesOutput = sprintf([str_subFolderRigid 'Vertices' str_Target str_BoneType str_BoneSide '.mat']);   
    [F,Y,~] = fImportSTL(str_TargetSTL);
    fnew=removeisolatedsurf(Y,F,10);
    [Y,F]=removeisolatednode(Y,fnew);
    %     figure;plotsurf(Y,F);
    %     figure;plotsurf(Y,fnew);
    %     figure;plotsurf(no,el);
    if strcmp(str_BoneType,'Patella')
        str_subfolderTargetSTL = [str_folder SubjectNumber '\\stlNew\\TransformStaticforCutNew\\'];
        str_MovingSTL = sprintf([str_subfolderTargetSTL 'DynamicPatella' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
        [F_long,X_long,~] = fImportSTL(str_MovingSTL);
        X_Cut = X_long;
        F_Cut = F_long;
    else
        str_subfolderTargetSTL = [str_folder SubjectNumber '\\stlNew\\TransformStaticforCutNew\\'];
        str_subfolderRefSTL = [str_folder SubjectNumber '\\stlNew\\Transformnocut\\'];
        if strcmp(str_BoneType,'Tibia')
            str_MovingSTL = sprintf([str_subfolderTargetSTL 'DynamicTibia' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_long,X_long,~] = fImportSTL(str_MovingSTL);
            str_RefSTL = sprintf([str_subfolderRefSTL 'DynamicFemur' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_Ref,X_Ref,~] = fImportSTL(str_RefSTL);
        else
            str_MovingSTL = sprintf([str_subfolderTargetSTL 'DynamicFemur' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_long,X_long,~] = fImportSTL(str_MovingSTL);
            str_RefSTL = sprintf([str_subfolderRefSTL 'DynamicTibia' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_Ref,X_Ref,~] = fImportSTL(str_RefSTL);
        end
        [X_Cut,F_Cut] = get_cut_bone(X_long,F_long, Y,F,X_Ref,lower(str_BoneType));
    end
    
    %% rigid transformation
    
    figure(1); clf();
    
    % Set the options
    opt.method='rigid';  % use rigid registration
    opt.corresp=1;       % estimate the correspondence vector
    opt.normalize=0;     % normalize to unit variance and zero mean before registering (default)
    opt.max_it=1000;      % max number of iterations
    opt.tol=1e-10;        % tolerance
    opt.viz=1;           % show every iteration
    opt.outliers=outliers;    % use 0.5 noise weight
    opt.fgt=fgt;           % [0,1,2] if > 0, then use FGT. case 1: FGT with fixing sigma after it gets too small (faster, but the result can be rough)
    %  case 2: FGT, followed by truncated Gaussian approximation (can be quite slow after switching to the truncated kernels, but more accurate than case 1)
    
    % Rigid registration options
    opt.rot=1;           % 1 - estimate strictly rotation. 0 - also allow for reflections.
    opt.scale=0;         % 1- estimate scaling. 0 - fixed scaling.
    
    % registering X to Y
    [Rigid, C]=cpd_register(Y,X_Cut,opt);
    
    figure(2); clf(); cpd_plot_iter(X_Cut, Y); title('Before');
    figure; cpd_plot_iter(Y, Rigid.Y); title('After rigid tranformation');
    
    % [Rigid, ~]=cpd_register(Y,X,opt);
    %
    % figure(2); clf(); cpd_plot_iter(Y, X); title('Before');
    % figure(3); cpd_plot_iter(Rigid.Y, Y); title('After rigid tranformation');
    
    RigidS.R = Rigid.R;
    RigidS.t = Rigid.t;
    RigidS.s = Rigid.s;
    save(str_RigidFileOutput,'RigidS')
    [dis_vertices,dist_surf,dist_vertices] = registration_error(Y,F,Rigid.Y,F_Cut, C);
    save(str_DistanceVOutput,'dis_vertices');
    save(str_DistancesurfOutput,'dist_surf');
    save(str_DistanceCorrespondenceOutput,'C');    
    save(str_DistanceVerticesOutput,'dist_vertices');   
    
    [regstr_matrix] = convertCPD2Affine3D(RigidS);
    fv.vertices = transformPointsForward(regstr_matrix,X_Cut);
    fv.faces = F_Cut;
    str_Output = sprintf([str_folder SubjectNumber '\\stlNew\\TransformPreNew\\Dynamic' str_BoneType str_BoneSide '_%d.stl'],nFrameNum-1);
    stlwrite(str_Output,fv);
    
    fv.vertices = transformPointsForward(regstr_matrix,X_long);
    fv.faces = F_long;
    str_Output = sprintf([str_folder SubjectNumber '\\stlNew\\TransformStaticforCutNew\\Dynamic' str_BoneType str_BoneSide '_%d.stl'],nFrameNum-1);
    stlwrite(str_Output,fv);
end

for nFrameNum = 5:10%7:11%
    str_Target = sprintf('D%dtoD%d',nFrameNum,nFrameNum+1);
    str_subFolderRigid = [str_folder SubjectNumber '\\matlab\\TransformPreNew\\'];
    str_RigidFileOutput = sprintf([str_subFolderRigid 'rigid' str_Target str_BoneType str_BoneSide '.mat']);
    
    % str_subfolderTargetSTL = [str_folder SubjectNumber '\\stl\\Static\\'];
    % str_MovingSTL = sprintf([str_subfolderTargetSTL str_BoneType '_' str_BoneSide '_reduce_cut.stl']);
    
    
    figure;plotsurf(X_long,F_long);
    % str_subfolderMovingSTL = [str_folder SubjectNumber '\\stl\\Dynamic\\'];
    str_TargetSTL = sprintf([str_DynamicFolder SubjectNumber '\\stlNew\\DynamicResample\\' str_BoneType str_BoneSide '_'  num2str(nFrameNum+1) '.stl']);
    str_DistanceVOutput = sprintf([str_subFolderRigid 'Dist_V' str_Target str_BoneType str_BoneSide '.mat']);
    str_DistancesurfOutput = sprintf([str_subFolderRigid 'Dist_surf' str_Target str_BoneType str_BoneSide '.mat']);
    str_DistanceCorrespondenceOutput = sprintf([str_subFolderRigid 'Correspondence' str_Target str_BoneType str_BoneSide '.mat']);    
    str_DistanceVerticesOutput = sprintf([str_subFolderRigid 'Vertices' str_Target str_BoneType str_BoneSide '.mat']);   
    [F,Y,~] = fImportSTL(str_TargetSTL); %%Hamid_FEMUR_Reduce.stl
    fnew=removeisolatedsurf(Y,F,10);
    [Y,F]=removeisolatednode(Y,fnew);
    %     figure;plotsurf(Y,F);
    %     figure;plotsurf(Y,fnew);
    %     figure;plotsurf(no,el);
    if strcmp(str_BoneType,'Patella')
        str_subfolderTargetSTL = [str_folder SubjectNumber '\\stlNew\\TransformStaticforCutNew\\'];
        str_MovingSTL = sprintf([str_subfolderTargetSTL 'DynamicPatella' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
        [F_long,X_long,~] = fImportSTL(str_MovingSTL);
        X_Cut = X_long;
        F_Cut = F_long;
    else
        str_subfolderTargetSTL = [str_folder SubjectNumber '\\stlNew\\TransformStaticforCutNew\\'];
        str_subfolderRefSTL = [str_folder SubjectNumber '\\stlNew\\Transformnocut\\'];
        if strcmp(str_BoneType,'Tibia')
            str_MovingSTL = sprintf([str_subfolderTargetSTL 'DynamicTibia' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_long,X_long,~] = fImportSTL(str_MovingSTL);
            str_RefSTL = sprintf([str_subfolderRefSTL 'DynamicFemur' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_Ref,X_Ref,~] = fImportSTL(str_RefSTL);
        else
            str_MovingSTL = sprintf([str_subfolderTargetSTL 'DynamicFemur' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_long,X_long,~] = fImportSTL(str_MovingSTL);
            str_RefSTL = sprintf([str_subfolderRefSTL 'DynamicTibia' str_BoneSide '_'  num2str(nFrameNum) '.stl']);
            [F_Ref,X_Ref,~] = fImportSTL(str_RefSTL);
        end
        [X_Cut,F_Cut] = get_cut_bone(X_long,F_long, Y,F,X_Ref,lower(str_BoneType));
    end
    %% rigid transformation
    
    figure(1); clf();
    
    % Set the options
    opt.method='rigid';  % use rigid registration
    opt.corresp=1;       % estimate the correspondence vector
    opt.normalize=0;     % normalize to unit variance and zero mean before registering (default)
    opt.max_it=1000;      % max number of iterations
    opt.tol=1e-10;        % tolerance
    opt.viz=1;           % show every iteration
    opt.outliers=outliers;    % use 0.5 noise weight
    opt.fgt=fgt;           % [0,1,2] if > 0, then use FGT. case 1: FGT with fixing sigma after it gets too small (faster, but the result can be rough)
    %  case 2: FGT, followed by truncated Gaussian approximation (can be quite slow after switching to the truncated kernels, but more accurate than case 1)
    
    % Rigid registration options
    opt.rot=1;           % 1 - estimate strictly rotation. 0 - also allow for reflections.
    opt.scale=0;         % 1- estimate scaling. 0 - fixed scaling.
    
    % registering X to Y
    [Rigid, C]=cpd_register(Y,X_Cut,opt);
    
    figure(2); clf(); cpd_plot_iter(X_Cut, Y); title('Before');
    figure; cpd_plot_iter(Y, Rigid.Y); title('After rigid tranformation');
    
    % [Rigid, ~]=cpd_register(Y,X,opt);
    %
    % figure(2); clf(); cpd_plot_iter(Y, X); title('Before');
    % figure(3); cpd_plot_iter(Rigid.Y, Y); title('After rigid tranformation');
    
    RigidS.R = Rigid.R;
    RigidS.t = Rigid.t;
    RigidS.s = Rigid.s;
    save(str_RigidFileOutput,'RigidS')
    [dis_vertices,dist_surf,dist_vertices] = registration_error(Y,F,Rigid.Y,F_Cut, C);
    save(str_DistanceVOutput,'dis_vertices');
    save(str_DistancesurfOutput,'dist_surf');
    save(str_DistanceCorrespondenceOutput,'C');    
     save(str_DistanceVerticesOutput,'dist_vertices'); 
    [regstr_matrix] = convertCPD2Affine3D(RigidS);
    fv.vertices = transformPointsForward(regstr_matrix,X_Cut);
    fv.faces = F_Cut;
    str_Output = sprintf([str_folder SubjectNumber '\\stlNew\\TransformPreNew\\Dynamic' str_BoneType str_BoneSide '_%d.stl'],nFrameNum+1);
    stlwrite(str_Output,fv);
    
    fv.vertices = transformPointsForward(regstr_matrix,X_long);
    fv.faces = F_long;
    str_Output = sprintf([str_folder SubjectNumber '\\stlNew\\TransformStaticforCutNew\\Dynamic' str_BoneType str_BoneSide '_%d.stl'],nFrameNum+1);
    stlwrite(str_Output,fv);
end

