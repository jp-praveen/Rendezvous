% Minimum time Rendezvous for a given Delta V
% Debri used: MICROSAT-R DEB 
% Debri TLE data: 1 44117U 19006C   19257.18743374  .00356745  20455-4  44493-2 0  9990
%                 2 44117  95.2886 104.5636 0572671 189.5932 169.3966 14.64658678 24390
% SGP4 output:  X                Y                Z     [km]
%            542.13411088    -3934.00754565    -5496.24138485 
%             XDOT             YDOT             ZDOT    [km/s]
%            -2.04410472        5.83405751       -4.77434023

setparameters(398588.738,104.5636,95.2886,189.5932,169.3966,0.572671)  % setting the global variable in the order (u,RAAN,Inclination,Perigee,Mean anomaly,e)
mat=getparameters;
u=mat(1,1);
e=mat(1,6);
RAAN=mat(1,2);
inclination=mat(1,3);
perigee=mat(1,4);

dv_max=8;

pos_t_ini=[542.13411088 -3934.00754565 -5496.24138485];          % in km
vel_t_ini=[-2.04410472 5.83405751 -4.77434023];                  % in kmpersec

ma_ini=mat(1,5);  %()
ta_ini=ma_ta(ma_ini);                                            % True anomaly of target's initial position

rn_t_ini=norm(pos_t_ini);
vn_t_ini=norm(vel_t_ini);
vr_t_ini=dot(vel_t_ini,pos_t_ini)/rn_t_ini;
vp_t_ini=sqrt(vn_t_ini^2-vr_t_ini^2);

h_t=rn_t_ini*vp_t_ini;                                  %specific angluar momentum of chaser

a_target=(h_t^2/u)/(1-e^2);                                    % in km
a_chaser=a_target+50;                                          % in km  
h_c=sqrt(a_chaser*u*(1-e^2));                           %specific angluar momentum of chaser

rp_target=a_target*(1-e);                               % Periapsis distance in target orbit 
rp_chaser=a_chaser*(1-e);                               % Periapsis distance in target orbit

T_t=(2*pi*a_target^(3/2))/u^0.5;
T_c=(2*pi*a_chaser^(3/2))/u^0.5;
i=-1;
for twait=0:100:T_c;                                   % Waiting Time
    i=i+1;
    ma_chaser=(2*pi*twait)/T_c;                        % mean anomaly of chaser in radians
    ma_chaser_deg=ma_chaser*180/pi;

    ta_chaser_deg=ma_ta(ma_chaser_deg);                % true anomaly of chaser in degrees   

    r1=tanomaly_rad(ta_chaser_deg,h_c);
    
    for dt=50:100:T_t;
        t=twait+dt;
        ma_target=(2*pi*t)/T_t;                         % mean anomaly of target in radians
        ma_target_deg=ma_target*180/pi;

        ta_target_deg=ma_ta(ma_target_deg);             % true anomaly of  target in degrees    
        r2=tanomaly_rad(ta_target_deg,h_t);
        
        [v1_pro,v2_pro,RAAN_pro,inclination_pro,perigee_pro,ta_1_pro,ta_2_pro,v1_retro,v2_retro,RAAN_retro,inclination_retro,perigee_retro,ta_1_retro,ta_2_retro]=lambert(r1,r2,dt);
        r1n=norm(r1);
        r2n=norm(r2);
        
        DCM=[cosd(RAAN),-sind(RAAN),0;sind(RAAN),cosd(RAAN),0;0,0,1]*[1 0 0;0 cosd(inclination) -sind(inclination);0 sind(inclination) cosd(inclination)]*[cosd(perigee+ta_target_deg) -sind(perigee+ta_target_deg) 0; sind(perigee+ta_target_deg) cosd(perigee+ta_target_deg) 0; 0 0 1];
        vp_t=h_t/r2n;                                  % Perpendicular velocity of the target
        vr_t=(u/h_t)*e*sind(ta_target_deg);            % Radial Velocity of the target 
        v_target=[vr_t;vp_t;0];                                   % velocity vector of the target 
        v_target_eci=DCM*v_target;
        
        DCM=[cosd(RAAN),-sind(RAAN),0;sind(RAAN),cosd(RAAN),0;0,0,1]*[1 0 0;0 cosd(inclination) -sind(inclination);0 sind(inclination) cosd(inclination)]*[cosd(perigee+ta_chaser_deg) -sind(perigee+ta_chaser_deg) 0; sind(perigee+ta_chaser_deg) cosd(perigee+ta_chaser_deg) 0; 0 0 1];
        vp_c=h_c/r1n;                                  % Perpendicular velocity of the chaser
        vr_c=(u/h_c)*e*sind(ta_chaser_deg);            % Radial Velocity of the chaser 
        v_chaser=[vr_c;vp_c;0];                        % Velocity Vector of the target
        v_chaser_eci=DCM*v_chaser;
        
        
        dv1_pro=v1_pro-v_chaser_eci;
        dv2_pro=v2_pro-v_target_eci;
        dv_prograde_1=norm(dv1_pro)+norm(dv2_pro);
        
        dv1_retro=v1_retro-v_chaser_eci;
        dv2_retro=v2_retro-v_target_eci;
        dv_retrograde_1=norm(dv1_retro)+norm(dv2_retro);
            
        
            if dv_prograde_1<dv_max | dv_retrograde_1<dv_max;
                if dv_prograde_1<dv_retrograde_1;
                    dv(1,i+1)=dv_prograde_1;
                else 
                    dv(1,i+1)=dv_retrograde_1;
                end
                tmin(1,i+1)=t;
                tcoast(1,i+1)=twait;
                transfer_time(1,i+1)=dt;
                true_anomaly_target(1,i+1)=ta_target_deg;
                true_anomaly_chaser(1,i+1)=ta_chaser_deg;
                r_cha(:,i+1)=r1;
                r_tar(:,i+1)=r2;
                r_chaser_norm(:,i+1)=r1n;
                r_target_norm(:,i+1)=r2n;
                break;
            end
    end    
end 
[a,b]=size(tcoast);
for k=1:b;
    if tcoast(1,1)==0 & tmin(1,1)==0 ;
        tcoast(tcoast==0)=[];
        tmin(tmin==0)=[];
        transfer_time(transfer_time==0)=[];
    end
end
plot(tcoast,tmin,'b-o');
title('Coasting Period VS Tmin');
xlabel('Coasting time (in seconds)');
ylabel('Tmin [Coasting time+Transfer time] (in seconds)');  
