% Minimum time Rendezvous for a given Delta V
% Debri used: MICROSAT-R DEB 
% Debri TLE data: 1 44117U 19006C   19257.18743374  .00356745  20455-4  44493-2 0  9990
%                 2 44117  95.2886 104.5636 0572671 189.5932 169.3966 14.64658678 24390
% SGP4 output:  X                Y                Z     [km]
%            542.13411088    -3934.00754565    -5496.24138485 
%             XDOT             YDOT             ZDOT    [km/s]
%            -2.04410472        5.83405751       -4.77434023
n=10;
m=1;
method=3;             % 1--> Howard ; 2--> Prussing ; 3--> Lancaster
u=398588.738;                             % in km^3*s^-2 
%dv_max=0.5;                                 % in km/sec
re=0;
fname = 'data_check.txt';

% Open the TLE file and read TLE elements
fid = fopen(fname, 'r');

for i=1:2*n;
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
    a(1,i) = re+( u/(no(1,i)*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    rNo(1,i) = str2num(tline(64:68));                % Revolution Number at Epoch
end
fclose(fid);

% Input 'r' and 'v' values. The first row in position and velocity data
% corresponds to the chaser. Except for debri_position,t_initial and
% velocity_position all arrays are (1,n).

% ----------------------- WATCH OUT HERE --------------------------------
[x,y,z,xdot,ydot,zdot]=test_sgp4;
% ----------------------- PRONE TO ERROR -------------------------------
for i=1:n;
    if i==1; 
        r_chaser_1=[re+x(i,1);re+y(i,1);re+z(i,1)];                % Initial position of chaser
        v_chaser_1=[xdot(i,1);ydot(i,1);zdot(i,1)];       % Velocity of chaser at initial position
        h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));               % Specific angluar momentum
        T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);              % Time Period   
        true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));
        t_initial_chaser=(mean_anomaly(1,i)*T(1,i))/(2*pi);
    else
        debri_position_1(:,i-1)=[re+x(i,1);re+y(i,1);re+z(i,1)];          % Initial position of debris; nth column implies nth debri data
        debri_velocity_1(:,i-1)=[xdot(i,1);ydot(i,1);zdot(i,1)]; % Initial velocity of debris; nth column implies nth debri data
        h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));                      % specific angluar momentum   
        T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);                  % Time Period
        t_initial(1,i-1)=(mean_anomaly(1,i)*T(1,i))/(2*180);
        true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));
    end  
end
j=0;
for twait=50000:1:150000;
    twait
    i=1;
    T(i+1);
    
    j=j+1;
    %twait=twait+5;
    for twait=twait;                                   % Waiting Time
        
        r1_ini=r_chaser_1;
        v1_ini=v_chaser_1;
        [r1,v1,alpha,universal_anomaly_chaser]=find_r2_v2(r1_ini,v1_ini,twait,T(1));
        ma_deg_chaser=(u^2/h(1,1)^3)*(1-e(1,1)^2)^(1.5)*(t_initial_chaser+twait)*180/pi;
        ta_chaser=ma_ta(ma_deg_chaser,e(1,1));
        k=0;
        r1=transpose(r1);
        v1=transpose(v1);
        
        for dt=500:50:10000;  %
            
            k=k+1;
            t=twait+dt;
            r2_ini=debri_position_1(:,i);
            v2_ini=debri_velocity_1(:,i);
            [r2,v2,alpha,universal_anomaly_target]=find_r2_v2(r2_ini,v2_ini,t,T(i+1));
            check_r2(i,k)=norm(r2);
            check_ua_target(i,k)=universal_anomaly_target;
            if isnan(r2)==0 ;
                dt;             
                ma_deg_target=(u^2/h(1,i+1)^3)*(1-e(1,i+1)^2)^(1.5)*(t_initial(1,i)+t)*180/pi;
                ta_target=ma_ta(ma_deg_target,e(1,i+1));
                dtheta=-1*(ta_chaser-ta_target);
                %if dtheta<0;
                 %   dtheta=360+dtheta
                %end
                %function [dv_prograde_1,dv_retrograde_1]=main_code_book_sub1(r1,r2,dt,dtheta,v1,v2);
                if method==1;
                    [v1_prograde,v2_prograde,RAAN_prograde,inclination_prograde,perigee_prograde,true_anomaly_1_prograde,true_anomaly_2_prograde,v1_retrograde,v2_retrograde,RAAN_retrograde,inclination_retrograde,perigee_retrograde,true_anomaly_1_retrograde,true_anomaly_2_retrograde] = lambert_book(r1,r2,dt);
                elseif method==2;
                    [v1_prograde,v2_prograde,at_short,RAAN_short,inclination_short,perigee_short,v1_retrograde,v2_retrograde,at_long,RAAN_long,inclination_long,perigee_long] = lambert_prussing(r1,r2,dt,dtheta,m);
                else
                    r2=transpose(r2);
                    v2=transpose(v2);
                    %m=1;
                    dt=dt/86400;
                    [V1, V2, extremal_distances, exitflag] = lambert(r1, r2, dt, m, u);
                    v1_prograde=V1;
                    v2_prograde=V2;
                    v1_retrograde=V1;
                    v2_retrograde=V2;
                    dt=dt*86400;
                end
            
                dv1_pro=v1_prograde-v1;
                dv2_pro=v2-v2_prograde;
                dv_prograde_1=norm(dv1_pro)+norm(dv2_pro);
        
                dv1_retro=v1_retrograde-v1;
                dv2_retro=v2-v2_retrograde;
                dv_retrograde_1=norm(dv1_retro)+norm(dv2_retro); 
                %end
        
                %if dv_prograde_1<dv_max || dv_retrograde_1<dv_max;
                
                    if isnan(dv_prograde_1)==0;
                        dv(j,k)=dv_prograde_1;
                    elseif isnan(dv_retrograde_1)==0; 
                        dv(j,k)=dv_retrograde_1;
                    else
                        dv(j,k)=nan;
                    end
                    %a_short(i,k)=at_short;
                    %a_long(i,k)=at_long;
                    
                    if isnan(dv(j,k));
                        tmin(j,k)=nan;
                        tcoast(j,k)=nan;
                        transfer_time(j,k)=nan;
                        ta_debri(j,k)=nan;
                        ma_debri(j,k)=nan;
                    else
                        tmin(j,k)=t;
                        tcoast(j,k)=twait;
                        transfer_time(j,k)=dt;
                        ta_debri(j,k)=ta_target;
                        ma_debri(j,k)=ma_deg_target;
                    end
                        
                    %delta_theta(i,k)=dtheta;
                    %ma_debri(1,i)=ma_deg;
                    %break;
                %end
            else
                dv(j,k)=nan;
                tmin(j,k)=nan;
                tcoast(j,k)=nan;
                transfer_time(j,k)=nan;
                %dt
            end
        end
    end
    
end
        
        
%[a,b]=size(tcoast);
%for k=1:n-1;
 %   if tmin(k,1)==0 ;
  %      tcoast(tcoast(k,:)==0)=[];
   %     tmin(tmin(k,:)==0)=[];
        %transfer_time(transfer_time==0)=[];
    %end
%end
n1=size(dv);
x=[];
y=[];
z=[];
for i=1:n1(1);
    dv_plot=dv(i,:);
    dv_plot=dv_plot(~isnan(dv_plot));
    dt=transfer_time(i,:);
    dt=dt(~isnan(dt));
    coasting_time=tcoast(i,:);
    coasting_time=coasting_time(~isnan(coasting_time));
    coasting_time(coasting_time==0)=nan;
    dt(dt==0)=nan;
    x=padconcatenation(x,dt,2);
    y=padconcatenation(y,coasting_time,2);
    z=padconcatenation(z,dv_plot,2); 
    %plot3(x,y,z,'.');
    %figure(i);
    %plot(dt,dv_plot);  
    dt=[];
    coasting_time=[];
    dv_plot=[];
end
v_x=linspace(min(x),max(x),100);
v_y=linspace(min(y),max(y),100);
[xx,yy]=meshgrid(v_x,v_y);
zz=griddata(x,y,z,xx,yy);
mesh(xx,yy,zz)
shading flat;
hold on
title('Tranfer time Vs Coasting time Vs dv');
xlabel('Transfer time (in seconds)');
ylabel('Coasting time (in seconds)');
zlabel('dv (in km/s)');

for i=1:n1(1);
    temp=dv(i,:);
    mindv(i)=min(temp(temp>0));
    if mindv(i)==[]
        mindv(i)=nan;
    end
end

mindv

%for i=1:n-1;
 %   if mindv(i)<1.2;
  %      useful_dv(i)=mindv(i);
   % end
%end

