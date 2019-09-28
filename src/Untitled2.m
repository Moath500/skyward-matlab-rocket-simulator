dataFile = zeros(length(data_ascent.integration.t)+...
    length(data_para.integration(1).t)+length(data_para.integration(2).t),4);
[~, ~, Press,~]= atmosisa([data_ascent.interp.alt'; data_para.interp(1).alt';
    data_para.interp(2).alt']);

%% ASCENT
% --  Acceleration in inertial frame
t = data_ascent.state.T;
N = length(t);
ua = data_ascent.state.Y(:,4);
% Compute acceleration with finite differences
aa = (ua(3:N)-ua(1:N-2))./(t(3:N)-t(1:N-2));
% Add derivatevs at boundaries
aa = [ua(2)/t(2); aa; (ua(end)-ua(end-1))/(t(end)-t(end-1))];
% Add static acceleration (in BODY frame)
as = quatrotate(data_ascent.state.Y(:,10:13),[0 0 -9.806]);
% Coriolis term in derivates. Recall: a_B = a_N - omega x V 
cor = cross(data_ascent.state.Y(:,7:9),data_ascent.state.Y(:,4:6));
% Recall: a_B = a_N - omega x V
aa = aa - cor(:,1) + as(:,1);

%% DROUGUE (INERTIAL = BODY frame)
t = data_para.state(1).T;
N = length(t);
ud = data_para.state(1).Y(:,4);
ad1 = (ud(3:N)-ud(1:N-2))./(t(3:N)-t(1:N-2));
ad1 = [ud(2)/t(2); ad1; (ud(end)-ud(end-1))/(t(end)-t(end-1))];
% Add static acceleration 
as = [-9.806 0 0]';
ad1 = ad1 + as(1);

t = data_para.state(2).T;
N = length(t);
ud = data_para.state(2).Y(:,4);
ad2 = (ud(3:N)-ud(1:N-2))./(t(3:N)-t(1:N-2));
ad2 = [ud(2)/t(2); ad2; (ud(end)-ud(end-1))/(t(end)-t(end-1))];
as = [-9.806 0 0]';
ad2 = ad2 + as(1);

ax = [aa; ad1; ad2]./9.80665;

%% WRITE TO FILE
dataFile(:,1) = T;
dataFile(:,2) = [-data_ascent.state.Y(:,3); -data_para.state(1).Y(:,3); ...
    -data_para.state(2).Y(:,3)];
dataFile(:,3) = Press';
dataFile(:,4) = ax;

fileID = fopen('para.txt','w');
formatSpec = '%4.2f %4.2f %6.2f %6.2f \r\n';
fprintf(fileID,'%4s %4s %6s %6s \n','t[s]','h[m]','P[Pa]','ax[g]');
fprintf(fileID,formatSpec,dataFile');
fclose(fileID);