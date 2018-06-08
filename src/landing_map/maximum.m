function [p_M] = maximum(x,y)
p = [x; y];
p = sqrt(sum(p.^2, 1)); %norm of each vector
[~, it_M] = max(p);
p_M = [x(it_M) y(it_M)];
end