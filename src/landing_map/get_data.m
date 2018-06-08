function [x,y]=getdata(origin)
open(origin);
h = gcf;

axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children'); 

objTypes = get(dataObjs, 'Type'); 

x = get(dataObjs, 'YData');  
y = get(dataObjs, 'XData');

close(h);

end
