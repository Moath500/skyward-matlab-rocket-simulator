function [h1,h2,h3,h4,h5] = front(mid_col,max_col,x,y,col)

% [x, y] = parsing(x, y, 'front');

%mean
[p_m, ad] = mid(x, y);

% %max
[p_M] = maximum(x,y);

% find min and max angle of landing points
theta = atan2(y,x)*180/pi;
theta(theta<0) = 360-abs(theta(theta<0)); % convert 0 to 360
theta_m = min(theta);
theta_M = max(theta);

%Plot
Rm = norm(p_m);
RM = norm(p_M); 

hold on
h1=annulus(0, 0, 0, RM, max_col, deg2rad(theta_m-10), deg2rad(theta_M+10));
h2=annulus(0, 0, Rm-ad, Rm+ad, mid_col, deg2rad(theta_m-5), deg2rad(theta_M+5));
alpha(h1,0.5);
alpha(h2,0.5);
h3=plot(p_M(1), p_M(2),'o','MarkerSize',13,'Color',col,'Linewidth',1.5);
% plot(p_m(1),p_m(2),'k.','Markersize',21,'Color',[1 1 0])

h4=plot(x, y, '.','Color','k','MarkerSize',10);
% h4=plot(x, y, '.','Color',col,'MarkerSize',10);
h5=plot(0,0,'r.','MarkerSize',21);


grid on
hold off
axis equal
end