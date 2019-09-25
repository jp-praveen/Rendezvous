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

dv_max=18;

pos_t_ini=[542.13411088 -3934.00754565 -5496.24138485];          % in km
vel_t_ini=[-2.04410472 5.83405751 -4.77434023];                  % in kmpersec

ma_ini=mat(1,5);  %()
ta_ini=ma_ta(ma_ini);                                            % True anomaly of target's initial position

rn_t_ini=norm(pos_t_ini);
vn_t_ini=norm(vel_t_ini);
vr_t_ini=dot(vel_t_ini,pos_t_ini)/rn_t_ini;
vp_t_ini=sqrt(vn_t_ini^2-vr_t_ini^2);
h_t=rn_t_ini*vp_t_ini;

a_target=(h_t^2/u)/(1-e^2);                                    % in km
a_chaser=a_target+50;                                          % in km  
h_c=sqrt(a_chaser*u*(1-e^2));                        %()

rp_target=a_target*(1-e);                               % Periapsis distance in target orbit 
rp_chaser=a_chaser*(1-e);                               % Periapsis distance in target orbit

T_t=(2*pi*a_target^(3/2))/u^0.5;
T_c=(2*pi*a_chaser^(3/2))/u^0.5;
i=-1;
for twait=0.1:100:0.25*T_c;                                   % Waiting Time
    i=i+1;
    ma_chaser=(2*pi*twait)/T_c;                        % mean anomaly of chaser in radians
    ma_chaser_deg=ma_chaser*180/pi;
    ta_chaser_deg=ma_ta(ma_chaser_deg);
    r1=tanomaly_rad(ta_chaser_deg,h_c);
            
    for dt=50:100:T_t;
        t=twait+dt;
        ma_target=(2*pi*t)/T_t;                         % mean anomaly of target in radians
        ma_target_deg=ma_target*180/pi;
        ta_target_deg=ma_ta(ma_target_deg);
        r2=tanomaly_rad(ta_target_deg,h_t);
        
        r1n=norm(r1);
        r2n=norm(r2);
        c=r2-r1;
        cn=norm(c);
        s1=(r1n+r2n+cn);
        s2=(r1n+r2n-cn);
        sets1s2(s1,s2,u,dt);
        
        a0=(a_target+a_chaser);
        at=fzero(@solve_at,a0);                            % Semi major axis of the transfer orbit
        
        if imag(at)==0 & at>s2/4;
            alpha_rad=2*asin(sqrt(s1/(4*at)));                   
            beta_rad=2*asin(sqrt(s2/(4*at)));
            alpha=alpha_rad*180/3.14;
            beta=beta_rad*180/3.14;
            g1=(1-r1n/at)*cotd(alpha-beta)-(1-r2n/at)*cscd(alpha-beta);
            g2=1-r1n/at;
        
            Eanomaly_trans_1=atand(g1/g2);
            Eanomaly_trans_2=Eanomaly_trans_1+alpha-beta;
            e_t=(1-r1n/at)/cosd(Eanomaly_trans_1);
            Tanomaly_trans_1=2*atand(sqrt((1+e_t)/(1-e_t))*tand(Eanomaly_trans_1/2));
            Tanomaly_trans_2=2*atand(sqrt((1+e_t)/(1-e_t))*tand(Eanomaly_trans_2/2));
            [RAAN_trans,inclination_trans,perigee_trans]=orbitparameters(r1,r2,e_t,Eanomaly_trans_1,Eanomaly_trans_2);
            DCM_trans=[cosd(RAAN_trans),-sind(RAAN_trans),0;sind(RAAN_trans),cosd(RAAN_trans),0;0,0,1]*[1 0 0;0 cosd(inclination_trans) -sind(inclination_trans);0 sind(inclination_trans) cosd(inclination_trans)]*[cosd(perigee_trans) -sind(perigee_trans) 0; sind(perigee_trans) cosd(perigee_trans) 0; 0 0 1];
        
            h_trans= sqrt(r1n*u*(1+e_t*cosd(Tanomaly_trans_1)));
            h_chaser=sqrt(u*r1n*(1+e));
            vp1_transfer=h_trans/r1n;                                               % Transfer orbit velocity at r1
            vr1_transfer=(u/h_trans)*e_t*sind(Tanomaly_trans_1);
            vp2_transfer=h_trans/r2n;                                               % Transfer orbit velocity at r1
            vr2_transfer=(u/h_trans)*e_t*sind(Tanomaly_trans_2);
            v1_transfer=[vr1_transfer;vp1_transfer;0];
            v1_transfer_eci=DCM_trans*v1_transfer;
            v2_transfer=[vr2_transfer;vp2_transfer;0];
            v2_transfer_eci=DCM_trans*v2_transfer;
        
            DCM=[cosd(RAAN),-sind(RAAN),0;sind(RAAN),cosd(RAAN),0;0,0,1]*[1 0 0;0 cosd(inclination) -sind(inclination);0 sind(inclination) cosd(inclination)]*[cosd(perigee) -sind(perigee) 0; sind(perigee) cosd(perigee) 0; 0 0 1];
        
            vp_t=h_t/r2n;                                  % Perpendicular velocity of the target
            vr_t=(u/h_t)*e*sind(ta_target_deg);            % Radial Velocity of the target 
            vt=[vr_t;vp_t;0];                        % velocity vector of the target
            vt_eci=DCM*vt;
        
            vp_c=h_c/r1n;                                  % Perpendicular velocity of the chaser
            vr_c=(u/h_c)*e*sind(ta_chaser_deg);            % Radial Velocity of the chaser 
            vc=[vr_c;vp_c;0];                        % Velocity Vector of the chaser
            vc_eci=DCM*vc;
        
            dv1=v1_transfer_eci-vc_eci;
            dv2=v2_transfer_eci-vt_eci;
            dv=norm(dv1)+norm(dv2);
        
            if dv<dv_max;
                deltav(1,i+1)=dv;
                tmin(1,i+1)=t;
                tcoast(1,i+1)=twait;
                transfer_time(1,i+1)=dt;
                true_anomaly_target(1,i+1)=ta_target_deg;
                r_chaser(:,i+1)=r1;
                r_target(:,i+1)=r2;
                trans_velocity_1(:,i+1)=v1_transfer_eci;
                trans_velocity_2(:,i+1)=v2_transfer_eci;
                true_anomaly_chaser(1,i+1)=ta_chaser_deg;
                
                break;
            end
        end
    end    
end 
tmin(tmin==0)=[];
transfer_time(transfer_time==0)=[];
tcoast(tcoast==0)=[];
deltav(deltav==0)=[];
true_anomaly_target(true_anomaly_target==0)=[];
true_anomaly_chaser(true_anomaly_chaser==0)=[];
r_chaser(r_chaser==0)=[];
r_target(r_target==0)=[];
trans_velocity_1(trans_velocity_1==0)=[];
trans_velocity_2(trans_velocity_2==0)=[];
plot(tcoast,tmin,'b-o');
title('Coasting Period VS Tmin');
xlabel('Coasting time (in seconds)');
ylabel('Tmin [Coasting time+Transfer time] (in seconds)');  
