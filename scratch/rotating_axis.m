function rotating_axis(R,limits)



% limits: vector containing 6 boundaries for axis and as last term the time
% for representation

if nargin==1
    x1=  -1;
    x2=   1;
    y1=  -1;
    y2=   1;
    z1=  -1;
    z2=   1;
    time=0.01;
else
    x1=limits(1);
    x2=limits(2);
    y1=limits(3);
    y2=limits(4);
    z1=limits(5);
    z2=limits(6);
    time=limits(7);
end

if size(R,3)>1
    [R]=square_to_line(R);
end

u1=R(:,1);
u2=R(:,2);
u3=R(:,3);
v1=R(:,4);
v2=R(:,5);
v3=R(:,6);
w1=R(:,7);
w2=R(:,8);
w3=R(:,9);

figure

for i=1:length(u1)
    
    clf('reset')
    grid on
    axis on
    axis([x1 x2 y1 y2 z1 z2])
    hold on
    quiver3(0,0,0,u1(i),v1(i),w1(i))
    quiver3(0,0,0,u2(i),v2(i),w2(i))
	quiver3(0,0,0,u3(i),v3(i),w3(i))
    hold off
    
    legend('x','y','z')

    
    pause(time)
    
end

end

function[R_l]=square_to_line(R_s)

n=size(R_s,3);
R_l=zeros(n,9);
for i=1:n
   
    R_l(i,:)=[R_s(1,1,i),R_s(2,1,i),R_s(3,1,i),R_s(1,2,i),R_s(2,2,i),R_s(3,2,i),R_s(1,3,i),R_s(2,3,i),R_s(3,3,i)];
    
end


end