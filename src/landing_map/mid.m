function [p_m, ad] = mid(x, y)
    
x_m = mean(x); 
y_m = mean(y);
p_m = [x_m y_m]; %mean point

p_n = [x; y];
p_n = sqrt(sum(p_n.^2, 1)); %norm2 of each vector(x,y)
ad = mad(p_n); %mean absolute deviation

end