function h = annulus(cx, cy, r1, r2, col, theta0, theta_end)
%PLOTANNULUS
%
%   H = PLOTANNULUS(CX, CY, R1, R2, COL, [THETA0])
%
%   THETA0 Start angle in radians of the first color, if more than one color
%       for the annulus is given. Default is 0.
%
%  Examples
%     plotannulus(0, 0, 5, 10, repmat([0 0 1; 1 1 0], 2, 1), deg2rad(45))
%     plotannulus(0, 0, 5, 10, parula(12))
%     plotannulus(0, 0, 5, 10, jet)
if nargin < 5, col = [0.5 0.5 0.5]; end
if nargin < 6, theta0 = 0; theta_end = 2*pi; end
Theta = linspace(theta0, theta_end);


  x1 = r1*cos(Theta);
  y1 = r1*sin(Theta);
    x2 = r2*cos(Theta);
    y2 = r2*sin(Theta);
    h = patch([x1 fliplr(x2)] + cx,...
                 [y1 fliplr(y2)] + cy,...
                 col, 'EdgeColor', col);
if nargout == 0
  clear h
end