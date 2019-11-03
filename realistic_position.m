function [time_elapsed,dv,next_node] = realistic_position(totalnode,startnode,twait,sequence_best)        
tn=totalnode;
s=startnode;
method=2;
u=398588.738;                             % in km^3*s^-2 
                                 
fname = 'data_check_paper_1to50.txt';

% Open the TLE file and read TLE elements
fid = fopen(fname, 'r');

for i=1:2*tn;
    % read first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    Cnum(1,i) = str2num(tline(3:7));      			        % Catalog Number (NORAD)
    epoch(1,i) = str2num(tline(19:32));              % Epoch
    Etype(1,i) = str2num(tline(63));                          % Ephemeris Type
    Enum(1,i)  = str2num(tline(65:end));             % Element Number
    
    % read second line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    inclination(1,i) = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    RAAN(1,i) = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e(1,i) = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    perigee(1,i) = str2num(tline(35:42));              % Argument of Perigee (degrees)
    mean_anomaly(1,i) = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    no(1,i) = str2num(tline(53:63));                 % Mean Motion
    a(1,i) =( u/(no(1,i)*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    rNo(1,i) = str2num(tline(64:68));                % Revolution Number at Epoch
end
fclose(fid);

[x,y,z,xdot,ydot,zdot]=test_sgp4;

for i=1:tn;
        debri_position_1(:,i)=[x(i,1);y(i,1);z(i,1)];          % Initial position of debris; nth column implies nth debri data
        debri_velocity_1(:,i)=[xdot(i,1);ydot(i,1);zdot(i,1)]; % Initial velocity of debris; nth column implies nth debri data
        h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));                      % specific angluar momentum   
        T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);                  % Time Period
        t_initial(1,i)=(mean_anomaly(1,i)*T(1,i))/(2*180);
        true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));
     
end

for i=1:tn
    r=debri_position_1(:,i);
    v=debri_velocity_1(:,i);
    [r1,v1,alpha,universal_anomaly]=find_r2_v2(r,v,twait,T(1,i));
    debri_position_new(i,:)=r1;
    debri_velocity_new(i,:)=v1;
end

r1=debri_position_new(s,:);
v1=debri_velocity_new(s,:);
r1=transpose(r1);
v1=transpose(v1);
ma_deg_chaser=(u^2/h(s)^3)*(1-e(s)^2)^(1.5)*(t_initial(s)+twait)*180/pi;
ta_chaser=ma_ta(ma_deg_chaser,e(s));

for i=1:tn;
    i
    if i~=sequence_best;
        r2_ini=debri_position_new(i,:);
        v2_ini=debri_velocity_new(i,:);
        k=0;
        for dt=100:100:2*T(i);
            k=k+1;
            ma_deg_target=(u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5)*(t_initial(1,i)+twait+dt)*180/pi;
            ta_target=ma_ta(ma_deg_target,e(1,i));
            dtheta=abs(ta_chaser-ta_target);
            [r2,v2,alpha,universal_anomaly_target]=find_r2_v2(r2_ini,v2_ini,dt,T(i));
            r2=transpose(r2);
            v2=transpose(v2);
            if method==1;
                [v1_prograde,v2_prograde,RAAN_prograde,inclination_prograde,perigee_prograde,true_anomaly_1_prograde,true_anomaly_2_prograde,v1_retrograde,v2_retrograde,RAAN_retrograde,inclination_retrograde,perigee_retrograde,true_anomaly_1_retrograde,true_anomaly_2_retrograde] = lambert_book(r1,r2,dt);
            else
                [v1_prograde,v2_prograde,at_short,RAAN_short,inclination_short,perigee_short,v1_retrograde,v2_retrograde,at_long,RAAN_long,inclination_long,perigee_long] = lambert_prussing(r1,r2,dt,dtheta);
            end
            
            dv1_pro=v1_prograde-v1;
            dv2_pro=v2-v2_prograde;
            dv_prograde_1=norm(dv1_pro)+norm(dv2_pro);
        
            dv1_retro=v1_retrograde-v1;
            dv2_retro=v2-v2_retrograde;
            dv_retrograde_1=norm(dv1_retro)+norm(dv2_retro); 
            
            if isnan(dv_prograde_1)==0;
                dv(i,k)=dv_prograde_1;
                transfer_time(i,k)=dt;
            elseif isnan(dv_retrograde_1)==0; 
                dv(i,k)=dv_retrograde_1;
                transfer_time(i,k)=dt;
            else
                dv(i,k)=nan;
                transfer_time(i,k)=nan;
            end
        end
    else 
        dv(i,1)=1000;
        transfer_time(i,1)=1000;
    end
end
%dv
%transfer_time
for i=1:tn;
    if i~=s;
        temp=dv(i,:);
        [minimum,index]=min(temp(temp>0));
        mindv(i)=minimum;
        mindv_dt(i)=transfer_time(i,index);
    else
        mindv(i)=100;
        mindv_dt(i)=100;
    end
end
dv=[];
transfer_time=[];
[dv,index]=min(mindv(mindv>0));
transfer_time=mindv_dt(index);
time_elapsed=transfer_time+twait;
next_node=index;
    
    






             
    
    