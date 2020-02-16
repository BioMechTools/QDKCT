 function [x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog, ... 
          a, R0, R] = lscylinder(X, x0, a0, r0, tolp, tolg)
% LSCYLINDER.M   Least-squares cylinder using Gauss-Newton.
%
% Version 1.0    
% Last amended   I M Smith 27 May 2002. 
% Created        I M Smith 08 Mar 2002
% Minor revisions  2016 MaxB
% ---------------------------------------------------------------------
% Input    
% X        Array [x y z] where x = vector of x-coordinates, 
%          y = vector of y-coordinates and z = vector of z-coordinates.
%          Dimension: m x 3. 
% 
% x0       Estimate of the point on the axis. 
%          Dimension: 3 x 1. 
%
% a0       Estimate of the axis direction. 
%          Dimension: 3 x 1. 
% 
% r0       Estimate of the cylinder radius. 
%          Dimension: 1 x 1. 
% 
% tolp     Tolerance for test on step length. 
%          Dimension: 1 x 1. 
%
% tolg     Tolerance for test on gradient.
%          Dimension: 1 x 1. 
% 
% Output  
% x0n      Estimate of the point on the axis. 
%          Dimension: 3 x 1. 
% 
% an       Estimate of the axis direction. 
%          Dimension: 3 x 1. 
% 
% rn       Estimate of the cylinder radius. 
%          Dimension: 1 x 1. 
% 
% d        Vector of radial distances from the points
%          to the cylinder. 
%          Dimension: m x 1. 
% 
% sigmah   Estimate of the standard deviation of the weighted 
%          residual errors. 
%          Dimension: 1 x 1. 
% 
% conv     If conv = 1 the algorithm has converged, 
%          if conv = 0 the algorithm has not converged
%          and x0n, rn, d, and sigmah are current estimates. 
%          Dimension: 1 x 1. 
% 
% Vx0n     Covariance matrix of point on the axis. 
%          Dimension: 3 x 3. 
%
% Van      Covariance matrix of axis direction. 
%          Dimension: 3 x 3. 
%
% urn      Uncertainty in cylinder radius. 
%          Dimension: 1 x 1. 
% 
% GNlog    Log of the Gauss-Newton iterations. 
%          Rows 1 to niter contain 
%          [iter, norm(f_iter), |step_iter|, |gradient_iter|]. 
%          Row (niter + 1) contains 
%          [conv, norm(d), 0, 0]. 
%          Dimension: (niter + 1) x 4. 
% 
% a        Optimisation parameters at the solution.
%          Dimension: 5 x 1. 
% 
% R0       Fixed rotation matrix. 
%          Dimension: 3 x 3. 
% 
% R        Upper-triangular factor of the Jacobian matrix
%          at the solution. 
%          Dimension: 5 x 5. 
%
% Modular structure: NLSS11.M, GNCC2.M, FGCYLINDER.M, ROT3Z.M, GR.M, 
%                    FGRROT3.M, FRROT3.M, DRROT3.M. 
%
% [x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog, a, R0, R] = 
%   lscylinder(X, x0, a0, r0, tolp, tolg <, w >)
% ---------------------------------------------------------------------


% check number of data points 
  m = size(X, 1);
  if m < 5
    error('At least 5 data points required: ' )
  end
% 
% find the centroid of the data 
  xb = mean(X)'; 
% 
% transform the data to close to standard position via a rotation 
% followed by a translation 
  R0 = rot3z(a0); % U * a0 = [0 0 1]' 
  x1 = R0 * x0; 
  xb1 = R0 * xb; 
% find xp, the point on axis nearest the centroid of the rotated data 
  t = x1 + (xb1(3) - x1(3)) * [0 0 1]'; 
  X2 = (X * R0') - ones(m ,1) * t'; 
  x2 = x1 - t; 
  xb2 = xb1 - t; 
% 
  ai = [0 0 0 0 r0]'; 
  tol = [tolp; tolg]'; 
% 
% Gauss-Newton algorithm to find estimate of roto-translation 
% parameters that transform the data so that the best-fit circle is 
% one in standard position
  [a, d, R, GNlog] = nlss11LOCAL(ai, tol, X2); 
% 
% inverse transformation to find axis and point on axis 
% corresponding to original data 
  rn = a(5); 
  [R3, DR1, DR2, DR3] = fgrrot3([a(3:4); 0]); 
  an = R0' * R3' * [0 0 1]'; % axis 
  p = R3 * (xb2 - [a(1) a(2) 0]'); 
  pz = [0 0 p(3)]'; 
  x0n = R0' * (t + [a(1) a(2) 0]' + R3' * pz); 
% x0n = point on axis in plane containing centroid of data 
% 
%   nGN = size(GNlog, 1); 
  conv = GNlog(end, 1); 
  if conv == 0 
    warning('*** Gauss-Newton algorithm has not converged ***'); 
  end % if conv 
% 
% calculate statistics 
  dof = m - 5; 
  sigmah = norm(d)/sqrt(dof); 
  ez = [0 0 1]'; 
  G = zeros(7, 5); 
% derivatives of x0n 
  dp1 = R3 * [-1 0 0]'; 
  dp2 = R3 * [0 -1 0]'; 
  dp3 = DR1 * (xb2 - [a(1) a(2) 0]'); 
  dp4 = DR2 * (xb2 - [a(1) a(2) 0]'); 
  G(1:3, 1) = R0' * ([1 0 0]' + R3' * [0 0 dp1'*ez]'); 
  G(1:3, 2) = R0' * ([0 1 0]' + R3' * [0 0 dp2'*ez]'); 
  G(1:3, 3) = R0' * (DR1' * [0 0 p'*ez]' + R3' * [0 0 dp3'*ez]'); 
  G(1:3, 4) = R0' * (DR2' * [0 0 p'*ez]' + R3' * [0 0 dp4'*ez]'); 
% derivatives of an 
  G(4:6, 3) = R0' * DR1' * [0 0 1]'; 
  G(4:6, 4) = R0' * DR2' * [0 0 1]'; 
% derivatives of rn 
  G(7, 5)   = 1; 
  Gt = R'\(sigmah * G'); % R' * Gt = sigmah * G' 
  Va = Gt' * Gt; 
  Vx0n = Va(1:3, 1:3); % covariance matrix for x0n 
  Van = Va(4:6, 4:6); % covariance matrix for an 
  urn = sqrt(Va(7, 7)); % uncertainty in rn 
% ---------------------------------------------------------------------
% End of LSCYLINDER.M.

function [a, f, R, GNlog] = nlss11LOCAL(ai,tol,X2)
%copy of nlss11, hardcoded fitting of cylinders to prevent using eval,
%'cause its evil. Also prevent reallocaiton of the GNlog
a0 = ai;
n = length(a0);
%
if n == 0
    error('Empty vector of parameter estimates:')
end
%
mxiter = (100+ceil(sqrt(n)));
conv = 0;
niter = 0;
eta = 0.01;
GNlog = zeros(mxiter+1,4) ; ii = 1 ;
%
% G-N iterations
while niter < mxiter && conv == 0
    [f0, J] = fgcylinder(a0,X2) ;%unit weigths are default
    if niter == 0 % scale by norm of columns of J
        nJ = size(J,2);
        scale = zeros(nJ,1);
        for j = 1:nJ;
            scale(j) = norm(J(:,j));
        end
    end % if niter
    %
    m = length(f0);
    % Check on m, n.
    if niter == 0 && m < n
        error('Number of observation less than number of parameters')
    end
    %
    % Calculate update step and gradient.
    F0 = norm(f0);
    Ra = triu(qr([J, f0]));
    R = Ra(1:nJ,1:nJ);
    q = Ra(1:nJ,nJ+1);
    p = -R\q;
    g = 2*R'*q;
    G0 = g'*p;
    a1 = a0 + p;
    niter = niter + 1;
    %
    % Check on convergence.
    [f1] = fgcylinder(a1,X2) ; %unit weights are default
    F1 = norm(f1);
    [c, conv, sp, sg] = gncc2(F0, F1, p, g, scale, tol(1), tol(2));
    %
    if conv ~= 1
        %
        % ...otherwise check on reduction of sum of squares.
        %
        % Evaluate f at a1.
        rho = (F1 - F0)*(F1 + F0)/(G0);
        if rho < eta
            tmin = max([0.001; 1/(2*(1-rho))]);
            a0 = a0 + tmin*p;
        else
            a0 = a0 + p;
        end % if rho
        %
    end % if conv
    %ii= 1
    GNlog(ii,:) = [niter, F0, sp, sg];%[GNlog; [niter, F0, sp, sg]];
    ii = ii + 1 ;
end % while niter
%
a = a0+p;
f = f1;
%
GNlog(ii,:) = [conv, F1, 0, 0];
GNlog = GNlog(1:ii,:) ;
%
% --------------------------------------------------------------------------
% End of NLSS11.M.

  function [c, iconv, sp, sg] = gncc2(f0, f1, p, g, scale, tolr, scalef)
% --------------------------------------------------------------------------
% GNCC2   Gauss-Newton convergence conditions.
%
% Version 1.0    
% Last amended   I M Smith 27 May 2002. 
% Created        I M Smith 08 Mar 2002
% --------------------------------------------------------------------------
% Input 
% f0       Norm of residuals at old point: norm(f(a)) 
%          Dimension: 1 x 1.
% 
% f1       Norm of residuals at new point: norm(f(a + p))
%          Dimension: 1 x 1.
%
% p        Step
%          Dimension: n x 1.
% 
% g        Gradient
%          Dimension: n x 1.
% 
% scale    Scale for columns of Jacobian matrix.
%          Dimension: n x 1.
% 
% tolr     Relative tolerance.
%          Dimension: 1 x 1.
% 
% scalef   Scale for function values. 
%          Dimension: 1 x 1.
%
% Output 
% c        Convergence indices in the form of ratio of value over
%          scaled tolerance.
%          c(1) size of step
%          c(2) change in function value
%          c(3) size of gradient
%          c(4) sum of squares near zero
%          c(5) gradient near zero
%          Dimension: 1 x 5.
% 
% iconv    = 1 if convergence criteria have been met, i.e., 
%          C(1), C(2), C(3) < 1 or 
%                      C(4) < 1 or 
%                      C(5) < 1.
%          = 0 otherwise.
%          Dimension: 1 x 1.
%
% sp       Scaled size of the step.
%          Dimension: 1 x 1.
% 
% sg       Scaled size of the gradient.
%          Dimension: 1 x 1.
%
% [c, iconv, sp, sg] = gncc2(f0, f1, p, g, scale, tolr, scalef)
% -------------------------------------------------------------------------- 

  iconv = 0;
%
  sp = max(abs(p .* scale));
  sg = max(abs(g ./ scale));
%
  c(1) = sp/(scalef * tolr^(0.7));
%
  delf = f0 - f1;
  c(2) = abs(delf)/(tolr * scalef);
%
  d3 = (tolr^(0.7)) * (scalef);
  d4 = scalef * (eps^(0.7));
  d5 = (eps^(0.7)) * (scalef);
%
  c(3) = sg/d3;
  c(4) = f1/d4;
  c(5) = sg/d5;
%
  if c(1) < 1 &  c(2) < 1 & c(3) < 1
    iconv = 1;
  elseif (c(4)) < 1 
    iconv = 1;
  elseif (c(5) < 1)
    iconv = 1;
  end
% --------------------------------------------------------------------------
% End of GNCC2.M.


function [f, J] = fgcylinder(a, X, w)
% ---------------------------------------------------------------------
% FGCYLINDER.M   Function and gradient calculation for
%                least-squares cylinder fit.
%
% Version 1.0
% Last amended   I M Smith 27 May 2002.
% Created        I M Smith 08 Mar 2002
% ---------------------------------------------------------------------
% Input
% a        Parameters [x0 y0 alpha beta s]'.
%          Dimension: 5 x 1.
%
% X        Array [x y z] where x = vector of x-coordinates,
%          y = vector of y-coordinates and z = vector of z-coordinates.
%          Dimension: m x 3.
%
% <Optional...
% w        Weights.
%          Dimension: m x 1.
% ...>
%
% Output
% f       Signed distances of points to cylinder:
%         f(i) = sqrt(xh(i)^2 + yh(i)^2) - s, where
%         [xh yh zh]' = Ry(beta) * Rx(alpha) * ([x y z]' - [x0 y0 0]').
%         Dimension: m x 1.
%
% <Optional...
% J       Jacobian matrix df(i)/da(j).
%         Dimension: m x 5.
% ...>
%
% Modular structure: FGRROT3.M, FRROT3.M, DRROT3.M.
%
% [f <, J >] = fgcylinder(a, X <, w >)
% ---------------------------------------------------------------------

m = size(X, 1);
% if no weights are specified, use unit weights
if nargin == 2
    w = ones(m, 1);
end % if nargin
%
x0 = a(1);
y0 = a(2);
alpha = a(3);
beta = a(4);
s = a(5);
%
[R, DR1, DR2] = fgrrot3([alpha beta 0]');
%
Xt = (X - ones(m, 1) * [x0 y0 0]) * R';
xt = Xt(:, 1);
yt = Xt(:, 2);
rt = sqrt(xt.*xt + yt.*yt);
Nt = [xt./rt yt./rt zeros(m, 1)];
f = dot(Xt, Nt, 2);
f = f - s;
f = w .* f; % incorporate weights
%
if nargout > 1 % form the Jacobian matrix
    J = zeros(m, 5);
    A1 = ones(m, 1) * (R * [-1 0 0]')';
    J(:, 1) = dot(A1, Nt, 2);
    A2 = ones(m, 1) * (R * [0 -1 0]')';
    J(:, 2) = dot(A2, Nt, 2);
    A3 = (X - ones(m, 1) * [x0 y0 0]) * DR1';
    J(:, 3) = dot(A3, Nt, 2);
    A4 = (X - ones(m, 1) * [x0 y0 0]) * DR2';
    J(:, 4) = dot(A4, Nt, 2);
    J(:, 5) = -1 * ones(m, 1);
    W = diag(w);
    J = W * J; % incorporate weights
end % if nargout
% ---------------------------------------------------------------------
% End of FGCYLINDER.M.

  function [U] = rot3z(a)
% --------------------------------------------------------------------------
% ROT3Z.M   Form rotation matrix U to rotate the vector a to a point along
%           the positive z-axis. 
%
% Version 1.0   
% Last amended  I M Smith 2 May 2002. 
% Created       I M Smith 2 May 2002. 
% --------------------------------------------------------------------------
% Input 
% a        Vector.
%          Dimension: 3 x 1. 
%
% Output 
% U        Rotation matrix with U * a = [0 0 z]', z > 0. 
%          Dimension: 3 x 3. 
% 
% Modular structure: GR.M. 
%
% [U] = rot3z(a)
% --------------------------------------------------------------------------

% form first Givens rotation
  [W, c1, s1] = gr(a(2), a(3));
  z = c1*a(2) + s1*a(3);
  V = [1 0 0; 0 s1 -c1; 0 c1 s1];
%
% form second Givens rotation
  [W, c2, s2] = gr(a(1), z);
%
% check positivity
  if c2 * a(1) + s2 * z < 0
    c2 = -c2;
    s2 = -s2;
  end  % if
%
  W = [s2 0 -c2; 0 1 0; c2 0 s2];
  U = W * V;
% --------------------------------------------------------------------------
% End of ROT3Z.M.

 function [U, c, s] = gr(x, y)
% --------------------------------------------------------------------------
% GR.M   Form Givens plane rotation. 
%
% Version 1.0    
% Last amended   I M Smith 27 May 2002. 
% Created        I M Smith 08 Mar 2002
% --------------------------------------------------------------------------
% Input 
% x        Scalar.
%          Dimension: 1 x 1. 
% 
% y        Scalar.
%          Dimension: 1 x 1. 
% 
% Output 
% U        Rotation matrix [c s; -s c], with U * [x y]' = [z 0]'. 
%          Dimension: 2 x 2. 
% 
% c        Cosine of the rotation angle. 
%          Dimension: 1 x 1. 
% 
% s        Sine of the rotation angle. 
%          Dimension: 1 x 1. 
% 
% [U, c, s] = gr(x, y)
% --------------------------------------------------------------------------

% form sine and cosine: s and c
  if y == 0
    c = 1;
    s = 0;
  elseif  abs(y) >= abs(x)
    t = x/y;
    s = 1/sqrt(1 + t*t);
    c = t*s;
  else
    t = y/x;
    c = 1/sqrt(1 + t*t);
    s = t*c;
  end % if
%
% assign U
  U = [c  s; -s  c];
% --------------------------------------------------------------------------
% End of GR.M.


  function [R, DR1, DR2, DR3] = fgrrot3(theta, R0)
% --------------------------------------------------------------------------
% FGRROT3.M   Form rotation matrix R = R3*R2*R1*R0 and its derivatives 
%             using right-handed rotation matrices:
%
%             R1 = [ 1  0   0 ]  R2 = [ c2 0  s2 ] and R3 = [ c3 -s3 0 ]
%                  [ 0 c1 -s1 ],      [ 0  1   0 ]          [ s3  c3 0 ].
%                  [ 0 s1  c2 ]       [-s2 0  c2 ]          [  0   0 1 ]
%
% Version 1.0    
% Last amended   I M Smith 27 May 2002. 
% Created        I M Smith 08 Mar 2002
% --------------------------------------------------------------------------
% Input 
% theta    Array of plane rotation angles (t1, t2, t3).
%          Dimension: 3 x 1. 
% 
% <Optional...  
% R0       Rotation matrix, optional, with default R0 = I.
%          Dimension: 3 x 3. 
% ...>
%
% Output 
% R        Rotation matrix.
%          Dimension: 3 x 3. 
% 
% <Optional...  
% DR1      Derivative of R wrt t1.
%          Dimension: 3 x 3. 
% 
% DR2      Derivative of R wrt t2.
%          Dimension: 3 x 3. 
% 
% DR3      Derivative of R wrt t3.
%          Dimension: 3 x 3. 
% ...>
%
% Modular structure: FRROT3.M, DRROT3.M. 
%
% [R <, DR1, DR2, DR3 >] = fgrrot3(theta <, R0 >)
% --------------------------------------------------------------------------

  if nargin == 1
    R0 = eye(3);
  end % if
%
  [R, R1, R2, R3] = frrot3(theta, R0);
%
% Evaluate the derivative matrices if required.
  if nargout > 1
	%used to b drrot3
    dR1 = [0 0 0; 0 -R1(3, 2) -R1(2, 2); 0 R1(2, 2) -R1(3, 2)];
    dR2 = [-R2(1, 3) 0 R2(1, 1); 0 0 0; -R2(1, 1) 0 -R2(1, 3)];
    dR3 = [-R3(2, 1) -R3(1, 1) 0; R3(1, 1) -R3(2, 1) 0; 0 0 0];
    DR1 = R3*R2*dR1*R0;
    DR2 = R3*dR2*R1*R0;
    DR3 = dR3*R2*R1*R0;
  end % if
% --------------------------------------------------------------------------
% End of FGRROT3.M.


  function [R, R1, R2, R3] = frrot3(theta, U0)
% --------------------------------------------------------------------------
% FRROT3.M   Form rotation matrix R = R3*R2*R1*U0. - use right-handed
%            rotation matrices.
%
% Version 1.0    
% Last amended   I M Smith 27 May 2002. 
% Created        I M Smith 08 Mar 2002
% --------------------------------------------------------------------------
% Input 
% theta    Array of plane rotation angles (t1, t2, t3).
%          Dimension: 3 x 1. 
% 
% <Optional...  
% U0       Rotation matrix, optional, with default R0 = I.
%          Dimension: 3 x 3. 
% ...>
%
% Output 
% R        Rotation matrix. 
%          Dimension: 3 x 3. 
% 
% R1       Plane rotation [1 0 0; 0 c1 -s1; 0 s1 c1].
%	       Dimension: 3 x 3. 
% 
% <Optional...  
% R2       Plane rotation [c2 0 s2; 0 1 0; -s2 0 c2].
%          Dimension: 3 x 3. 
% 
% R3       Plane rotation [c3 -s3 0; s3 c3 0; 0 0 1].
%          Dimension: 3 x 3. 
% ...>
%
% [R, R1 <, R2, R3 >] = frrot3(theta <, U0 >)
% --------------------------------------------------------------------------

  ct = cos(theta); 
  st = sin(theta);
%
  if length(theta) > 0
    R1 = [ 1 0 0; 0 ct(1) -st(1); 0 st(1) ct(1)];
    R = R1;
  end %if
%
  if length(theta) > 1
    R2 = [ ct(2) 0 st(2); 0 1 0; -st(2) 0 ct(2)];
    R = R2*R;
  end % if
%
  if length(theta) > 2
    R3 = [ ct(3) -st(3) 0; st(3) ct(3) 0; 0 0 1];
    R = R3*R;
   end % if
%
  if nargin > 1
    R = R*U0;
  end % if
% --------------------------------------------------------------------------
% End of FRROT3.M.