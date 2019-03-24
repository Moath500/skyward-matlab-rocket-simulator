function [h1,h2,h3,h4] = front(x,y)

%mean
[Rm, ad] = mid(x, y);

%max
[RM] = maximum(x,y);

% find min and max angle of landing points
theta = atan2(y,x)*180/pi;
theta(theta<0) = 360-abs(theta(theta<0)); % convert 0 to 360
theta_m = min(theta);
theta_M = max(theta);

hold on
mid_col = [0 128/255 1];
max_col = [204/255 229/255 1];
h1 = annulus(0, 0, 0, RM, max_col, deg2rad(theta_m), deg2rad(theta_M));
h2 = annulus(0, 0, Rm-ad, Rm+ad, mid_col, deg2rad(theta_m), deg2rad(theta_M));
alpha(h1,0.5);
alpha(h2,0.5);

h3 = plot(x, y, '.','Color','k','MarkerSize',10);
h4 = plot(0,0,'r.','MarkerSize',21);


grid on
hold off
axis equal
end

function [RM] = maximum(x,y)
N = length(x);
R = zeros(N,1);

for i = 1:N
R(i) = norm([x(i),y(i)]);
RM = max(R);
end

end

function [Rm, ad] = mid(x, y)
N = length(x);
R = zeros(N,1);

for i = 1:N
R(i) = norm([x(i),y(i)]);
Rm = mean(R);
end

ad = std(R);
end

function h = annulus(cx, cy, r1, r2, col, theta0, theta_end)

if nargin < 5, col = [0.5 0.5 0.5]; end
if nargin < 6, theta0 = 0; theta_end = 2*pi; end
Theta = linspace(theta0, theta_end);


  x1 = r1*cos(Theta);
  y1 = r1*sin(Theta);
    x2 = r2*cos(Theta);
    y2 = r2*sin(Theta);
    h = patch([x1 fliplr(x2)] + cx, [y1 fliplr(y2)] + cy, col, 'EdgeColor', col);
if nargout == 0
  clear h
end
end