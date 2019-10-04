% GIVEN INITAL POSITIONS AND TRANSFER TIME FIND THE VELOCITIES OF THE PROGRADE TRANSFER ORBIT
function [v1_prograde,v2_prograde,RAAN_prograde,inclination_prograde,perigee_prograde,true_anomaly_1_prograde,true_anomaly_2_prograde,v1_retrograde,v2_retrograde,RAAN_retrograde,inclination_retrograde,perigee_retrograde,true_anomaly_1_retrograde,true_anomaly_2_retrograde] = lambert(r1,r2,transfer_time)
u=398588.738;                             % in km^3*s^ 
dt=transfer_time;
r1cr2=cross(r1,r2);
r1r2=norm(r1)*norm(r2);
r1n=norm(r1);
r2n=norm(r2);

%PROGRADE ORBIT TRANSFER
if r1cr2(3,1)>=0;
    deltheta=acosd(dot(r1,r2)/r1r2);                                  % deltheta is change in true anomaly in transfer orbit or Transfer angle
else
    deltheta=360-acosd(dot(r1,r2)/r1r2);                              % deltheta is in degrees  
end    
A_prograde=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));
setsolveparameters(u,dt,r1n,r2n,A_prograde);
z0=10;

z_prograde=fsolve(@solve_torbit,z0);

if z_prograde>0 & imag(z_prograde)==0;
    S_prograde=(1/6)-(z_prograde/120)+(z_prograde^2/5040)-(z_prograde^3/362880)+(z_prograde^4/39916800);  % S and C are Stumpff Functions 
    C_prograde=(1/2)-(z_prograde/24)+(z_prograde^2/720)-(z_prograde^3/40320)+(z_prograde^4/3628800);
    y_prograde=r1n+r2n+A_prograde*((z_prograde*S_prograde-1)/C_prograde^0.5);

% Calculating the Lagrange functions and Velocities at r1 and r2

    f_prograde=1-y_prograde/r1n;
    g_prograde=(A_prograde*sqrt(y_prograde)/sqrt(u));
    gdot_prograde=1-(y_prograde/r2n);
    v1_prograde=(r2-f_prograde*r1)/g_prograde;                                               % Transfer orbit velocity at r1
    v2_prograde=(gdot_prograde*r2-r1)/g_prograde;                                           % Transfer orbit velocity at r2

% Calculating the Orbital Elements of the transfer orbit 
    evec = ((norm(v1_prograde)^2-u/norm(r1))*r1-dot(r1,v1_prograde)*v1_prograde)/u;
    e_prograde = norm(evec);

    v1_r=dot(r1,v1_prograde)/r1n;
    if v1_r>=0
        true_anomaly_1_prograde=acosd(dot(evec,r1)/(e_prograde*r1n));
    else
        true_anomaly_1_prograde=360-acosd(dot(evec,r1)/(e_prograde*r1n));
    end
    v2_r=dot(r2,v2_prograde)/r2n;
    if v2_r>=0
        true_anomaly_2_prograde=acosd(dot(evec,r2)/(e_prograde*r2n));
    else
        true_anomaly_2_prograde=360-acosd(dot(evec,r2)/(e_prograde*r2n));
    end
   
    h_prograde=cross(transpose(r1),transpose(v1_prograde));
    K=[0,0,1];
    n_prograde=cross(K,h_prograde);
    inclination_prograde=acosd((h_prograde(1,3)/norm(h_prograde)));
    RAAN_prograde = acosd(n_prograde(1,1)/norm(n_prograde));
    if n_prograde(1,2)<0;
        RAAN_prograde = 360-RAAN_prograde;
    end

    perigee_prograde = acosd(dot(n_prograde,evec)/(norm(n_prograde)*e_prograde));

    if evec(3,1)<0;
        perigee_prograde = 360-perigee_prograde;
    end
else
    v1_prograde=nan;
    v2_prograde=nan;
    RAAN_prograde=nan;
    inclination_prograde=nan;
    perigee_prograde=nan;
    true_anomaly_1_prograde=nan;
    true_anomaly_2_prograde=nan;
end


%RETROGRADE ORBIT TRANSFER

if r1cr2(3,1)>=0;
    deltheta=360-acosd(dot(r1,r2)/r1r2);                                  % deltheta is change in true anomaly in transfer orbit or Transfer angle
else
    deltheta=acosd(dot(r1,r2)/r1r2);                              % deltheta is in degrees  
end    
A_retrograde=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));
setsolveparameters(u,dt,r1n,r2n,A_retrograde);
z0=10;
 
z_retrograde=fsolve(@solve_torbit,z0);
if z_retrograde>0 & imag(z_retrograde)==0; 
    S_retrograde=(1/6)-(z_retrograde/120)+(z_retrograde^2/5040)-(z_retrograde^3/362880)+(z_retrograde^4/39916800);  % S and C are Stumpff Functions 
    C_retrograde=(1/2)-(z_retrograde/24)+(z_retrograde^2/720)-(z_retrograde^3/40320)+(z_retrograde^4/3628800);
    y_retrograde=r1n+r2n+A_retrograde*((z_retrograde*S_retrograde-1)/C_retrograde^0.5);
 
% Calculating the Lagrange functions and Velocities at r1 and r2
 
    f_retrograde=1-y_retrograde/r1n;
    g_retrograde=(A_retrograde*sqrt(y_retrograde)/sqrt(u));
    gdot_retrograde=1-(y_retrograde/r2n);
    v1_retrograde=(r2-f_retrograde*r1)/g_retrograde;                                               % Transfer orbit velocity at r1
    v2_retrograde=(gdot_retrograde*r2-r1)/g_retrograde;                                           % Transfer orbit velocity at r2
 
% Calculating the Orbital Elements of the transfer orbit
    evec = ((norm(v1_retrograde)^2-u/norm(r1))*r1-dot(r1,v1_retrograde)*v1_retrograde)/u;
    e_retrograde = norm(evec);
    v1_r=dot(r1,v1_retrograde)/r1n;
    if v1_r>=0
        true_anomaly_1_retrograde=acosd(dot(evec,r1)/(e_retrograde*r1n));
    else
        true_anomaly_1_retrograde=360-acosd(dot(evec,r1)/(e_retrograde*r1n));
    end
    v2_r=dot(r2,v2_retrograde)/r2n;
    if v2_r>=0
        true_anomaly_2_retrograde=acosd(dot(evec,r2)/(e_retrograde*r2n));
    else
        true_anomaly_2_retrograde=360-acosd(dot(evec,r2)/(e_retrograde*r2n));
    end
 
    h_retrograde=cross(transpose(r1),transpose(v1_retrograde));
    K=[0,0,1];
    n_retrograde=cross(K,h_retrograde);
    inclination_retrograde=acosd((h_retrograde(1,3)/norm(h_retrograde)));
    RAAN_retrograde = acosd(n_retrograde(1,1)/norm(n_retrograde)); 
    if n_retrograde(1,2)<0;
       RAAN_retrograde = 360-RAAN_retrograde;
    end
    perigee_retrograde = acosd(dot(n_retrograde,evec)/(norm(n_retrograde)*e_retrograde));
    if evec(3,1)<0;
        perigee_retrograde = 360-perigee_retrograde;
    end
else
    v1_retrograde=nan;
    v2_retrograde=nan;
    RAAN_retrograde=nan;
    inclination_retrograde=nan;
    perigee_retrograde=nan;
    true_anomaly_1_retrograde=nan;
    true_anomaly_2_retrograde=nan;
end
  


