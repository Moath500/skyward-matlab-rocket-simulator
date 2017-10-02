function [] = cone(origin,col)
%Parsing
[x, y, ~] = getdata(origin); 

[x, y] = parsing(x, y, 'front');

p_M = maximum(x,y);
R = norm(p_M);

[max_y, ind] = max(abs(y));
max_x = x(ind);
maxy = [max_x max_y];
x_proj = [max_x 0];

theta = acos(norm(x_proj)/norm(maxy));

hold on
annulus(0, 0, 0, R, col, -theta, theta)

plot(x, y, 'k*')
plot(0,0,'r.','MarkerSize',21)

grid on
hold off
axis equal
end