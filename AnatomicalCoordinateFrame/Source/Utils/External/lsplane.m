function [x0, a, d, normd] = lsplane(X)
% LSPLANE.M Least-squares plane (orthogonal distance regression).
%
% Input
% X Array [x y z] where x = vector of x-coordinates,
% y = vector of y-coordinates and z = vector of
% z-coordinates.
% Dimension: m x 3.
%
% Output
% x0 - Centroid of the data = point on the best-fit plane.
% Dimension: 3 x 1.
%
% a - Direction cosines of the normal to the best-fit plane.
% Dimension: 3 x 1.
%
% <Optional...
% d - Residuals.
% Dimension: m x 1.
%
% normd - Norm of residual errors.
% Dimension: 1 x 1.
% ...>
%
% [x0, a <, d, normd >] = lsplane(X)
%
% Created I M Smith 08 Mar 2002

% check number of data points
m = size(X, 1);
if m < 3
    error('At least 3 data points required: ' )
end

% calculate centroid
x0 = mean(X)';

% form matrix A of translated points
A = [(X(:, 1) - x0(1)) (X(:, 2) - x0(2)) (X(:, 3) - x0(3))];

% calculate the singular value decomposition of A
[U, S, V] = svd(A, 0);

% find the smallest singular value in S and extract from V the
% corresponding right singular vector
[s, i] = min(diag(S));
a = V(:, i);

% calculate residual distances, if required
if nargout > 2
    d = U(:, i)*s;
    normd = norm(d);
end