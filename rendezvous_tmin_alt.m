% Given initial set of debri orbit, find the least time taken by the chaser 
% to visit all the debri 
function [visit_sequence, min_t] = rendezvous_tmin_alt(n,x,y,z,xdot,ydot,zdot)

u=398588.738;                             % in km^3*s^-2 
dv_max=20;                                 % in km/sec
flag=[1:n-1];
% TLE file name 
fname = 'debri_data_swap.txt';

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
    a(1,i) = ( u/(no(1,i)*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    rNo(1,i) = str2num(tline(64:68));                % Revolution Number at Epoch
end
fclose(fid);

% Input 'r' and 'v' values. The first row in position and velocity data
% corresponds to the chaser. Except for debri_position,t_initial and
% velocity_position all arrays are (1,n).
for i=1:n;
    if i==1; 
        r_chaser_1=[x(i,1);y(i,1);z(i,1)];                % Initial position of chaser
        v_chaser_1=[xdot(i,1);ydot(i,1);zdot(i,1)];       % Velocity of chaser at initial position
        h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));               % Specific angluar momentum
        T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);              % Time Period   
        true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));
    else
        debri_position_1(:,i-1)=[x(i,1);y(i,1);z(i,1)];          % Initial position of debris; nth column implies nth debri data
        debri_velocity_1(:,i-1)=[xdot(i,1);ydot(i,1);zdot(i,1)]; % Initial velocity of debris; nth column implies nth debri data
        h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));                      % specific angluar momentum   
        T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);                  % Time Period
        t_initial(1,i-1)=(mean_anomaly(1,i)*T(1,i))/(2*pi);
        true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));
    end  
end

ta_chaser=true_anomaly(1,1);
k=10;
j=0;
while k~=0;
    j=j+1;
    for i=1:n-1;
        DCM=[cosd(RAAN(1,1)),-sind(RAAN(1,1)),0;sind(RAAN(1,1)),cosd(RAAN(1,1)),0;0,0,1]*[1 0 0;0 cosd(inclination(1,1)) -sind(inclination(1,1));0 sind(inclination(1,1)) cosd(inclination(1,1))]*[cosd(perigee(1,1)+ta_chaser) -sind(perigee(1,1)+ta_chaser) 0; sind(perigee(1,1)+ta_chaser) cosd(perigee(1,1)+ta_chaser) 0; 0 0 1];
        r1=DCM*r_chaser_1;
        v1=DCM*v_chaser_1;
        for dt=100:200:4000;
            r2_ini=debri_position_1(:,i);
            v2_ini=debri_velocity_1(:,i);
            [r2,v2,alpha]=find_r2_v2(r2_ini,v2_ini,dt);
            
            ma_deg=(u^2/h(1,i+1)^3)*(1-e(1,i+1)^2)^(1.5)*(t_initial(1,i)+dt);
            ta=ma_ta(ma_deg,e(1,i+1));
            DCM=[cosd(RAAN(1,i+1)),-sind(RAAN(1,i+1)),0;sind(RAAN(1,i+1)),cosd(RAAN(1,i+1)),0;0,0,1]*[1 0 0;0 cosd(inclination(1,i+1)) -sind(inclination(1,i+1));0 sind(inclination(1,i+1)) cosd(inclination(1,i+1))]*[cosd(perigee(1,i+1)+ta) -sind(perigee(1,i+1)+ta) 0; sind(perigee(1,i+1)+ta) cosd(perigee(1,i+1)+ta) 0; 0 0 1];
            r2=DCM*r2_ini;
            v2=DCM*v2_ini;
            [v1_short,v2_short,RAAN_short,inclination_short,perigee_short,v1_long,v2_long,RAAN_long,inclination_long,perigee_long] = lambert_prussing(r1,r2,dt);
            dv1_short=v1_short-v1;
            dv2_short=v2_short-v2;
            dv_short_1=norm(dv1_short)+norm(dv2_short);
        
            dv1_long=v1_long-v1;
            dv2_long=v2_long-v2;
            dv_long_1=norm(dv1_long)+norm(dv2_long);   
            if dv_short_1<dv_max || dv_long_1<dv_max;
                if dv_short_1<dv_long_1;
                    dv(1,i)=dv_short_1;
                else 
                    dv(1,i)=dv_long_1;
                end
                
                tmin(1,i)=dt;
                ta_debri(1,i)=ta;
                ma_debri(1,i)=ma_deg;
                break;
            end
        end
    end
    check=tmin;
    min_t(1,j)=min(tmin(tmin>0));
    s_vec=find(tmin==min_t(1,j));
    [o,p]=size(s_vec);
    
    for i=1:p;
        select(1,i)=flag(s_vec(i));
    end
    s1=min(select);
    s=find(select==s1);
    select=[];
        
    
% Updating Chasers new position, velocity and orbital elements 
    ta_chaser=ta_debri(1,s);
    inclination(1,1)=inclination(1,s+1);
    RAAN(1,1)=RAAN(1,s+1);
    e(1,1)=e(1,s+1);
    a(1,1)=a(1,s+1);
    T(1,1)=T(1,s+1);
    perigee(1,1)=perigee(1,s+1);
    h(1,1)=h(1,s+1);
    r_chaser_1=debri_position_1(:,s);
    v_chaser_1=debri_velocity_1(:,s);

% Updating the debri new position and velocities and removing the visited
% debris
    a(:,s+1)=[];
    h(:,s+1)=[];
    T(:,s+1)=[]; 
    ma_debri(:,s)=[];
    debri_position_1(:,s)=[];
    debri_velocity_1(:,s)=[];
    t_initial(:,s)=[];
    inclination(:,s+1)=[];
    RAAN(:,s+1)=[];
    e(:,s+1)=[];
    perigee(:,s+1)=[];
    
    [l,k]=size(debri_position_1);
    for i=1:k;
        t_initial(1,i)=(ma_debri(1,i)*T(1,i+1))/(2*pi);
        dt=min_t(1,j);
        r1=debri_position_1(:,i);
        v1=debri_velocity_1(:,i);
        [r2,v2,alpha]=find_r2_v2(r1,v1,dt);
        debri_position_1(:,i)=r2(:,1);
        debri_velocity_1(:,i)=v2(:,1);
        
    end
    visit_sequence(1,j)=flag(1,s);
    flag(:,s)=[];
    n=n-1;
    tmin=[];
    ta_debri=[];
end

end


            
       
        
      
        
    
    

