% Test windgen

for i=1:50
    [u,v,w]=windgen(pi/2,pi/2,0,0,0,20);
    quiver(u,v);
    hold on
end
