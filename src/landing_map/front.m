function [h1,h2,h3,h4,h5] = front(mid_col,max_col,x,y,col)

[x, y] = parsing(x, y, 'front');

%mean
[p_m, ad] = mid(x, y);

% %max
[p_M] = maximum(x,y);

%Plot
Rm = norm(p_m);
RM = norm(p_M); 

hold on
h1=annulus(0, 0, 0, RM, max_col, deg2rad(-90), deg2rad(90));
h2=annulus(0, 0, Rm-ad, Rm+ad, mid_col, deg2rad(-90), deg2rad(90));
alpha(h1,0.5);
alpha(h2,0.5);
h3=plot(p_M(1), p_M(2),'o','MarkerSize',13,'Color',col,'Linewidth',1.5);
% plot(p_m(1),p_m(2),'k.','Markersize',21,'Color',[1 1 0])


h4=plot(x, y, '.','Color',col,'MarkerSize',10);
h5=plot(0,0,'r.','MarkerSize',21);


grid on
hold off
axis equal
end