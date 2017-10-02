function [x, y ,z] = getdata(origin)
%Getting data vectors from figure
f = open(origin);
h = gcf;


axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');

x = get(dataObjs, 'XData');
y = get(dataObjs, 'YData');
z = get(dataObjs, 'ZData');

close(f)
end