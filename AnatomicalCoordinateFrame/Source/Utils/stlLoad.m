function [output] = stlLoad(varargin)
%STLLOAD load triangular mesh binary stl files.
%   [[output]] = stlLoad(filenames,[filespath],[0])
%   returns a struct, or, when no outputarguments are requested, saves it
%   as a .mat file. Ensures no duplicate vertices in the loaded data.
%
%---Input
%   filenames   -   (Cell array of) string(s) with name(s) of the STL files
%   filespath   -   Path to the stl files. Only one path is allowed
%   1           -   Add a 1to the argument list to draw the loaded files
%---Output
%    stlStruc   -   A struct array where the ith element has the fields F
%    and V holding the faces and vertices of the ith file from filenames
%
%    When no output is requested the things are saved in mat files.
%
%  Example usage (load all stl files in the current folder) :
%    allStlFiles = dir('*.stl')
%    meshes = stlLoad({allStlFiles.name},'.') ;
%
%  see also: dir save 

showPlot = any(cellfun(@(x) isnumeric(x) && x==1,varargin)) ;
% Open STL-files
%rm all numeric input from the input list
varargin = varargin(~cellfun(@isnumeric,varargin)) ;

if nargin >= 1
    stlfiles = varargin{1} ;
    if nargin < 2
        dashLocs = strfind(stlfiles,filesep) ;
        if isempty(dashLocs)
            stlpath = '.\' ;
        else
            stlpath = stlfiles(1:dashLocs(end)) ;
            stlfiles = stlfiles(dashLocs(end)+1:end ) ;
        end
    else
        stlpath = varargin{2} ;
    end
end
if iscellstr(stlfiles) == 1  
    nstlfile = length(stlfiles);
    stlfilesC = stlfiles ;
else
    nstlfile = 1;
    stlfilesC{1,1} = stlfiles;
    clear stlfiles;
end

if stlpath(end)~= filesep ;
    stlpath = [stlpath filesep];
end

for k = nstlfile:-1:1 %reverse iterate to prevent reallocation
    [F, V, N, stltitle] = stlreadbin_nested(stlfilesC{1,k},stlpath);
    %remove duplicate vertices
    [V,~,ic] = unique(V,'rows') ;
    F = ic(F) ;
    saveName{1,k} = [stlfilesC{1,k}(1:end-4),'.mat'];
    stlfile = stlfilesC{1,k};
    if nargout > 0
        output(k).F = F ;
        output(k).V = V ;
        output(k).N = N ;
        output(k).stltitle = stltitle ;
        output(k).stlfile = stlfile ;
        output(k).saveName = saveName{1,k} ;
    else
        % Save STLfile in -mat
        %         save([stlpath,SaveName{1,k}],'F','V','N','C','stlfile','-mat')
        save([stlpath,saveName{1,k}],'F','V','N','stlfile','-mat')
    end
end

if showPlot
    % Show figure of the loaded stl-file
    
    axis off,daspect([1,1,1]), view(3)
    
    hold on
    colorpoints = colormap('HSV');
    colorpoints = colorpoints(round(linspace(1,64,nstlfile)),:);
    for k = 1:nstlfile
        % Display stlfiles
        patch('Faces',output(k).F,'Vertices',output(k).V,'FaceColor',colorpoints(k,:),'EdgeColor','None','FaceAlpha',0.5);
        material dull
    end
    light; lighting phong; camlight;
end


function [F, V, N, stltitle,C] = stlreadbin_nested(stlfile,stlpath)
% This function reads an STL file in binary format into vertex and face
% matrices v and f.
%
% @@Maybe wanna rework this to allow bigger stls
%       Mostly inlining
% USAGE: [F, V, N, stltitle,C] = stlreadbin(stlfile,stlpath);
%
% F contains the vertex lists defining each triangle face [n x 3].
% V contains the vertices for all triangles [3*n x 3].
% N contains the normals for each triangle face [n x 3].
% stltitle contains the title of the specified stl file [1 x 80].
% C is optional and contains color rgb data in 5 bits [n x 3].
%
% To see plot the 3D surface use:
%   patch('Faces',F,'Vertices',V,'FaceVertexCData',c);
% or
%   plot3(V(:,1),V(:,2),V(:,3),'.');

% Based on code originally written by:
% Francis Esmonde-White, May 2010
% Re-written by WJ Zevenbergen, May 2012 (ORL Nijmegen)
% MaxChange: changed C to last output to prevent calculations

use_color = (nargout>=5);

% Select or check the stl-file
if nargin < 1;
    [stlfile, stlpath] = uigetfile('*.stl','Select STL-file');
elseif nargin == 1;
    % stlpath = Current Path
    if ~exist(stlfile,'file')
        error(['File ','%s',' not found.'], stlfile);
    end
else
    if ~exist([stlpath,stlfile],'file')
        error(['File ','%s',' not found.'], stlfile);
    end
end

fid = fopen([stlpath,stlfile], 'r'); % Open the file, assumes STL Binary format.
if fid == -1
    error('File could not be opened, check name or path')
end

% Read stl-file
ftitle = fread(fid,80,'uchar=>schar'); % Read file title
stltitle = char(ftitle');
numFaces = fread(fid,1,'int32'); % Read number of Faces
if numFaces == 0
    warning('No data in STL file');
    return
end

fprintf('\nTitle: %s\n', stltitle);
fprintf('Filename: %s\n', stlfile);
fprintf('Num Faces: %d\n', numFaces);

tic % measure time reading stl-file

T = fread(fid,inf,'uint8=>uint8'); % read the remaining values
fclose(fid);

% Each facet is 50 bytes
%  - Three single precision values specifying the face normal vector
%  - Three single precision values specifying the first vertex (XYZ)
%  - Three single precision values specifying the second vertex (XYZ)
%  - Three single precision values specifying the third vertex (XYZ)
%  - Two color bytes (possibly zeroed)

% 3 dimensions x 4 bytes x 4 vertices = 48 bytes for triangle vertices
% 2 bytes = color (if color is specified)

trilist = 1:48;

ind = reshape(repmat(50*(0:(numFaces-1)),[48,1]),[1,48*numFaces])+repmat(trilist,[1,numFaces]);
Tri = reshape(typecast(T(ind),'single'),[3,4,numFaces]);

N = squeeze(Tri(:,1,:))';
N = double(N);

V = Tri(:,2:4,:);
V = reshape(V,[3,3*numFaces]);
V = double(V)';

F = reshape(1:3*numFaces,[3,numFaces])';

if use_color
    c0 = typecast(T(49:50),'uint16');
    if (bitget(c0(1),16)==1)
        trilist = 49:50;
        ind = reshape(repmat(50*(0:(numFaces-1)),[2,1]),[1,2*numFaces])+repmat(trilist,[1,numFaces]);
        c0 = reshape(typecast(T(ind),'uint16'),[1,numFaces]);
        
        r = bitshift(bitand(2^16-1, c0),-10);
        g = bitshift(bitand(2^11-1, c0),-5);
        b = bitand(2^6-1, c0);
        C = [r; g; b]';
    else
        C = zeros(numFaces,3);
    end
end
readingtime = toc;
fprintf('STL reading time: %8.2f seconds\n',readingtime)


