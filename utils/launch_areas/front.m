function [] = front(origin,mid_col,max_col)
%Parsing
[x, y, ~] = getdata(origin); 

[x, y] = parsing(x, y, 'front');

%mean
[p_m, ad] = mid(x, y);

% %max
[p_M] = maximum(x,y);

%Plot
Rm = norm(p_m);
RM = norm(p_M); 

hold on
annulus(0, 0, 0, RM, max_col, deg2rad(-90), deg2rad(90))
annulus(0, 0, Rm-ad, Rm+ad, mid_col, deg2rad(-90), deg2rad(90))

plot(p_M(1), p_M(2),'o','MarkerSize',21,'Color',[0 0 102/255],'Linewidth',1.5)
plot(p_m(1),p_m(2),'k.','Markersize',21,'Color',[1 1 0])

plot(x, y, 'k*')
plot(0,0,'r.','MarkerSize',21)

grid on
hold off
axis equal
end