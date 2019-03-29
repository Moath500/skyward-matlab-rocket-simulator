function [h1,h2,h3,h4] = front(x,y)

%mean
[Rm, ad] = mid(x, y);

%max
[RM] = maximum(x,y);

% find min and max angle of landing points
theta = atan2(y,x);
theta_m = min(theta);
theta_M = max(theta);

hold on
mid_col = [0 128/255 1];
max_col = [204/255 229/255 1];
h1 = annulus(0, RM, max_col, theta_m, theta_M);
h2 = annulus(Rm-ad, Rm+ad, mid_col, theta_m, theta_M);
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

function h = annulus(r1, r2, col, theta0, theta_end)

Theta = linspace(theta0, theta_end);


x1 = r1*cos(Theta);
y1 = r1*sin(Theta);
x2 = r2*cos(Theta);
y2 = r2*sin(Theta);
h = patch([x1 fliplr(x2)], [y1 fliplr(y2)], col, 'EdgeColor', col);
if nargout == 0
    clear h
end
end