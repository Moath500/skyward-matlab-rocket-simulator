load para500.mat
norm_wind=zeros(100,1);
norm_R=zeros(100,1);
for i=1:100
    wind=data_para{i}.wind(1).body_wind(1:3);  
    norm_wind(i) = norm(wind);
    
    
    norm_R(i)= norm(data_para{i}.state(2).Y(end,1:2)); 
end

figure
plot(norm_wind,norm_R,'o')

[norm_wind,i]=sort(norm_wind);
norm_R=norm_R(i);


p=polyfit(norm_wind,norm_R,1);
y=polyval(p,norm_wind);

hold on

plot(norm_wind,y)


