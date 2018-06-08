function [x1, y1] = parsing(x, y, type)

l = length(x);
index = [];
if strcmp(type,'back')
    for k = 1:l
    if x(k)>0
        index = [index k];
    end
    end
elseif strcmp(type,'front')
    for k = 1:l
    if x(k)<0
        index = [index k];
    end
    end
else
    error('Select a type: front or back')
end

x(index) = [];
y(index) = [];

x1 = x;
y1 = y;

end