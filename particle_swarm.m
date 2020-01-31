iterations=1000;
points=10000;
coasting_time_1=linspace(0,150000,points);
transfer_time_1=linspace(500,15000,points);
coasting_time_2=coasting_time_1;
transfer_time_2=transfer_time_1;
velocity_ct_1=10*rand(1,points);
velocity_tt_1=10*rand(1,points);
dv_1=1500*ones(1,points);
position_ct=[];
position_tt=[];
velocity_tt_2=[];
velocity_ct_2=[];


n=2;
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

for i1=1:iterations;
    i1
    
    i=1;
    j=1;
    k=0;
    
       for i1=1:points;  %
            k=k+1;
            
            twait=coasting_time_2(i1);
            dt=transfer_time_2(i1);
            
            r1_ini=r_chaser_1;
            v1_ini=v_chaser_1;
            [r1,v1,alpha,universal_anomaly_chaser]=find_r2_v2(r1_ini,v1_ini,twait,T(1));
            ma_deg_chaser=(u^2/h(1,1)^3)*(1-e(1,1)^2)^(1.5)*(t_initial_chaser+twait)*180/pi;
            ta_chaser=ma_ta(ma_deg_chaser,e(1,1));
           
            r1=transpose(r1);
            v1=transpose(v1);
            
            
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
            end
       
                    
            if isnan(dv_prograde_1)==0;
                dv_2(j,k)=dv_prograde_1;
            elseif isnan(dv_retrograde_1)==0; 
                dv_2(j,k)=dv_retrograde_1;
            else
                dv_2(j,k)=1500;
            end
       end
       
       for i1=1:points
           if dv_2(1,i1)<dv_1(1,i1)
               dv(1,i1)=dv_2(1,i1);
               position_ct(1,i1)=coasting_time_2(1,i1);
               position_tt(1,i1)=transfer_time_2(1,i1);
           else
               dv(1,i1)=dv_1(1,i1);
               position_ct(1,i1)=coasting_time_1(1,i1);
               position_tt(1,i1)=transfer_time_1(1,i1);
           end
       end
                
       
       [dv_global,mark]=min(dv);
       dv_global
       position_ct(1,mark)
       position_tt(1,mark)
       for i1=1:points
           velocity_ct_2(1,i1)=0.4*velocity_ct_1(1,i1)+0.5*rand(1)*(position_ct(1,i1)-coasting_time_2(1,i1))+0.5*rand(1)*(position_ct(1,mark)-coasting_time_2(1,i1));
           velocity_tt_2(1,i1)=0.4*velocity_tt_1(1,i1)+0.5*rand(1)*(position_tt(1,i1)-transfer_time_2(1,i1))+0.5*rand(1)*(position_tt(1,mark)-transfer_time_2(1,i1));
       end
       
       coasting_time_1=coasting_time_2;
       transfer_time_1=transfer_time_2;
       dv_1=dv_2;
       velocity_ct_1=velocity_ct_2;
       velocity_tt_1=velocity_tt_2;
       
       
       for i1=1:points
           coasting_time_2(1,i1)=coasting_time_2(1,i1)+velocity_ct_2(1,i1);
           transfer_time_2(1,i1)=transfer_time_2(1,i1)+velocity_tt_2(1,i1);
       end
end
           
dv_global
position_ct(1,mark)
position_tt(1,mark)
       
           
           
                
                    