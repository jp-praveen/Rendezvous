% Minimum time Rendezvous for a given Delta V
% Debri used: MICROSAT-R DEB 
% Debri TLE data: 1 44117U 19006C   19257.18743374  .00356745  20455-4  44493-2 0  9990
%                 2 44117  95.2886 104.5636 0572671 189.5932 169.3966 14.64658678 24390
% SGP4 output:  X                Y                Z     [km]
%            542.13411088    -3934.00754565    -5496.24138485 
%             XDOT             YDOT             ZDOT    [km/s]
%            -2.04410472        5.83405751       -4.77434023
n=6;
method=2;
u=398588.738;                             % in km^3*s^-2 
dv_max=0.5;                                 % in km/sec
re=0;
fname = 'debri_data.txt';

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
[x,y,z,xdot,ydot,zdot]=test_sgp4;

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

for i=1:n-1;
    i
    T(i+1)
    j=0;
    
    for twait=1;                                   % Waiting Time
        j=j+1;
        r1_ini=r_chaser_1;
        v1_ini=v_chaser_1;
        [r1,v1,alpha]=find_r2_v2(r1_ini,v1_ini,twait);
        ma_deg_chaser=(u^2/h(1,1)^3)*(1-e(1,1)^2)^(1.5)*(t_initial_chaser+twait)*180/pi;
        ta_chaser=ma_ta(ma_deg_chaser,e(1,1));
        k=0;
        for dt=1:50:T(i+1);   %
            k=k+1;
            t=twait+dt;
            r2_ini=debri_position_1(:,i);
            v2_ini=debri_velocity_1(:,i);
            [r2,v2,alpha]=find_r2_v2(r2_ini,v2_ini,t);
            if isnan(r2)==0 ;
                
                
                ma_deg_target=(u^2/h(1,i+1)^3)*(1-e(1,i+1)^2)^(1.5)*(t_initial(1,i)+t)*180/pi;
                ta_target=ma_ta(ma_deg_target,e(1,i+1));
                dtheta=abs(ta_chaser-ta_target);
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
        
                %if dv_prograde_1<dv_max || dv_retrograde_1<dv_max;
                
                    if isnan(dv_prograde_1)==0;
                        dv(i,k)=dv_prograde_1;
                    elseif isnan(dv_retrograde_1)==0; 
                        dv(i,k)=dv_retrograde_1;
                    else
                        dv(i,k)=nan;
                    end
                    a_short(i,k)=at_short;
                    a_long(i,k)=at_long;
                    tmin(i,k)=t;
                    tcoast(i,k)=twait;
                    transfer_time(i,k)=dt;
                    ta_debri(i,k)=ta_target;
                    ma_debri(i,k)=ma_deg_target;
                    delta_theta(i,k)=dtheta;
                    %ma_debri(1,i)=ma_deg;
                    %break;
                %end
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
for i=1:n-1;
    dv_plot=dv(i,:);
    dt=transfer_time(i,:);
    coasting_time=tcoast(i,:);
    coasting_time(coasting_time==0)=nan;
    dt(dt==0)=nan;
    figure(i);
    plot(dt,dv_plot);
    grid
    title('Tranfer time Vs dv');
    xlabel('Transfer time (in seconds)');
    ylabel('dv (in km/s)');  
    dt=[];
    coasting_time=[];
end

