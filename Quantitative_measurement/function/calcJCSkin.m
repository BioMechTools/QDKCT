function [ varargout ] = calcJCSkin( varargin )
%CALCJCSKIN Kinematics in Joint Coordinate System
%   JCSkin = CALCJCSKIN(A, B, SIDE) calculate clinical rotations and 
%   translation using z-x-y sequence of rotations according to Grood and
%   Suntay (1983), where A, B are roto-translation matrices of two moving
%   reference frames in a common global frame.
%   
%   [JCSrot, JCStrans] = CALCJCSKIN(...) returns two N-by-3 matrices of 
%   rotations and translations.
%
%   Input:
%   A, B  4-by-4 position matrices of frames A and B (static configuration)
%         4-by-4-by-N position matrices of frames A and B (N time frames)
%         3-by-4-by-N position matrices (if last row [0 0 0 1] is omitted)
%         
%                     Xu Xv Xw Xo
%         u v w o  =  Yu Yv Yw Yo
%         0 0 0 1     Zu Zv Zw Zo
%                     0  0  0  1
%   
%   SIDE  case-insensitive string for side specification ('left', 'l',
%         'right', 'r')
%   
%   Output:
%   JCSkin  N-by-6 matrix (columns 1-3: rotations, 4-6: translations)
%           Rotations are the Tait-Bryan angles from body A (e.g. femur)
%           to body B (e.g. tibia, patella), using the sequence z-x-y:
%           1) flexion is around body A's z-axis
%           2) adduction is around a floating axis (x')
%           3) external rotation is around body B's y-axis
%           Translations are the projection of the vector from A to B onto
%           the joint coordinate system:
%           4) anterior displacement of body B
%           5) distal displacement of body B
%           6) lateral displacement of body B
%   

narginchk(2, 3);
nargoutchk(0, 2);

% check and uniformize number of frames
if size(varargin{1}, 3) ~= size(varargin{2}, 3)
  % input matrices have different sizes (one has only one frame)
  [NFRAMES, VARN]= max(size(varargin{1}, 3), size(varargin{2}, 3));
  switch VARN
    case 1
      A= varargin{1};
      B= repmat(varargin{2}, 1, 1, NFRAMES);
    case 2
      A= repmat(varargin{1}, 1, 1, NFRAMES);
      B= varargin{1};
  end
else
  % input matrices have same number of frames
  NFRAMES= size(varargin{1}, 3);
  A= varargin{1};
  B= varargin{2};
end

% both input matrices have the same size (3-by-3-by-N)
if size(varargin{1}, 1) == 3
  if size(varargin{1}, 2) == 3
    A= [A repmat([0;0;0], 1, 1, NFRAMES);
      repmat([0 0 0 1], 1, 1, NFRAMES)];
    B= [B repmat([0;0;0], 1, 1, NFRAMES);
      repmat([0 0 0 1], 1, 1, NFRAMES)];
  else
    if size(varargin{1}, 2) == 4
      A= [A; repmat([0 0 0 1], 1, 1, NFRAMES)];
      B= [B; repmat([0 0 0 1], 1, 1, NFRAMES)];
    else
      error('Wrong size of input matrices!');
    end
  end
end

if nargin < 3
  warning('No side specification provided.')
  signstr= input('Enter ''L'', for left or ''R'', for right: ', 's');
else
  signstr= varargin{3};
end  

if sum(strcmpi(signstr, {'left' 'l'}))
  SIGN= -1;
else
  if sum(strcmpi(signstr, {'right' 'r'}))
    SIGN= 1;
  else
    error('Wrong side specification! Possible options: ''right'', ''left''');
  end
end

T= zeros(3, 3, NFRAMES);
for iFrame= 1:NFRAMES
  T(:, :, iFrame)= B(1:3, 1:3, iFrame)'*A(1:3, 1:3, iFrame);
end

% convert rotation matrix to Euler angles (sequence Z-X-Y)
JCSrot= SpinCalc('DCMtoEA312', T, 1e-12, 0);

% sign manipulation to obtain flexion, adduction, tibial external rotation
JCSrot= JCSrot*[[-1 0 0]; [0 SIGN*1 0]; [0 0 -SIGN*1]];

JCStrans= zeros(NFRAMES, 3);
for iFrame= 1:NFRAMES
  % solve the equation AH = b, where H and b are vectors, A matrix
  H= A(:, :, iFrame)\B(:, 4, iFrame);
  
  alpha= JCSrot(iFrame, 1);
  beta= 90 + SIGN*JCSrot(iFrame, 2);
  
  V= [
    [cosd(alpha) -sind(alpha) 0]; ...
    [-sind(alpha)*sind(beta) -cosd(alpha)*sind(beta) -cosd(beta)]; ...
    [0 0 SIGN] ...
    ];
  
  q= V*H(1:3);
  JCStrans(iFrame, :)= q';
end

if nargout <= 1
  varargout{1}= [JCSrot JCStrans];
else
  varargout{1}= JCSrot;
  varargout{2}= JCStrans;
end

end

