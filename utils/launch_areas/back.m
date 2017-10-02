function [] = back(origin,mid_col,max_col)
%Parsing
[x, y, ~] = getdata(origin); 

[x, y] = parsing(x, y, 'back');

%mean
[p_m, ad] = mid(x, y);

% %max
[p_M] = maximum(x,y);

%Plot
Rm = norm(p_m);
RM = norm(p_M); 

hold on
annulus(0, 0, 0, RM, max_col, deg2rad(90), deg2rad(270))
annulus(0, 0, Rm-ad, Rm+ad, mid_col, deg2rad(90), deg2rad(270))

plot(p_M(1), p_M(2),'o','MarkerSize',21,'Color',[102/255 102/255 0],'Linewidth',1.5)
plot(p_m(1),p_m(2),'k.','Markersize',21,'Color',[1 128/255 0])

plot(x, y, 'k*')
plot(0,0,'r.','MarkerSize',21)

grid on
hold off
axis equal
end