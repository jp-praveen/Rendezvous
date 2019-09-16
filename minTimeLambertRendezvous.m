% Single Revolution Lambert Problem
% Debri used: MICROSAT-R DEB 
% Debri TLE data: 1 44117U 19006C   19257.18743374  .00356745  20455-4  44493-2 0  9990
%                 2 44117  95.2886 104.5636 0572671 189.5932 169.3966 14.64658678 24390
% Inclination (degrees):95.2886
% RAAN (degrees):104.5636
% Eccentricity:0.572671
% Argument of Perigee [Degrees]:189.5932
% Mean Anomaly [Degrees]:169.3966
% SGP4 output:  X                Y                Z     [km]
%            542.13411088    -3934.00754565    -5496.24138485 
%             XDOT             YDOT             ZDOT    [km/s]
%            -2.04410472        5.83405751       -4.77434023

e=0.572671;
u=398588.738;                                                % in km^3*s^-2                
pos_d=[542.13411088 -3934.00754565 -5496.24138485];          % in km
vel_d=[-2.04410472 5.83405751 -4.77434023];                  % in kmpersec
r_d=norm(pos_d);
v_d=norm(vel_d);
vr_d=dot(vel_d,pos_d)/r_d;
vp_d=sqrt(v_d^2-vr_d^2);
h_d=r_d*vp_d;
a_debri=(h_d^2/u)/(1-e^2);                                    % in km
a_target=a_debri+50;                                          % in km  
rp_debri=a_debri*(1-e);
rp_target=a_target*(1-e);
r1_cylindrical=rp_target;

%Calculating time period(T) of debri and finding the position at t=T/4
T=(2*pi*a_debri^(3/2))/u^0.5;
t_quarter=T/4;
ma_quarter=(2*pi*t_quarter)/T;                     % ma=mean anomaly in radians
if ma_quarter<pi;                                  % ea=eccentric anomaly     
    eai_quarter=ma_quarter+e/2;                    % initiating ea
else
    eai_quarter=ma_quarter-e/2;
end
       
f=eai_quarter-e*sin(eai_quarter)-ma_quarter;
f_dash=1-e*cos(eai_quarter);
while abs((f/f_dash))>0.000001;                         % Newtons method to find eccentric anomaly 
    ea_quarter=eai_quarter-(f/f_dash);
    eai_quarter=ea_quarter;
    f=eai_quarter-e*sin(eai_quarter)-ma_quarter;
    f_dash=1-e*cos(eai_quarter); 
end

ta_quarter=2*atand(sqrt((1+e)/(1-e))*tand(ea_quarter/2))*180/3.14;    % ta=True Anomaly (in degrees) 
r_quarter=(h_d^2/u)*(1+e*cosd(ta_quarter))^(-1);
r2_cylindrical=r_quarter;

% Converting radial vectors to ECI frame
ta_debri=ta_quarter;
ta_chaser=0;                                             % True Anomaly of chaser
r1_cartesian=[r1_cylindrical*cosd(ta_chaser);r1_cylindrical*sind(ta_chaser);0];
r2_cartesian=[r2_cylindrical*cosd(ta_debri);r1_cylindrical*sind(ta_debri);0];

o=104.5636;                                             % RAAN
i=95.2886;                                              % Inclination
w=189.5932;                                             % Argument of Periapsis
DCMc=[cosd(o),-sind(o),0;sind(o),cosd(o),0;0,0,1]*[1 0 0;0 cosd(i) -sind(i);0 sind(i) cosd(i)]*[cosd(w) -sind(w) 0; sind(w) cosd(w) 0; 0 0 1];
DCMt=[cosd(o),-sind(o),0;sind(o),cosd(o),0;0,0,1]*[1 0 0;0 cosd(i) -sind(i);0 sind(i) cosd(i)]*[cosd(w+ta_debri) -sind(w+ta_debri) 0; sind(w+ta_debri) cosd(w+ta_debri) 0; 0 0 1];
r1=DCMc*r1_cartesian;                                    % Writing r1 and r2 in ECI frame
r2=DCMt*r2_cartesian;

% Assuming both target and chaser are at periapsis at t=0 and choosing a prograde trajectory(i<90 deg). Rendezvous occurs at t=T/4 of debri orbit
r1cr2=cross(r1,r2);
r1r2=norm(r1)*norm(r2);
if r1cr2(3,1)>=0;
    deltheta=acosd(dot(r1,r2)/r1r2);
else
    deltheta=360-acosd(dot(r1,r2)/r1r2);
end    

A=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));           

% Starting Newtons method to evaluate the z (z is  related to a and universal anomaly) 
z=0;                                                     % Initializing z
S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);  % S and C are Stumpff Functions 
C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
r1n=norm(r1);
r2n=norm(r2);
del_t=t_quarter;
y=r1n+r2n+A*((z*S-1)/C^0.5);
f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
if z==0;
    f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
else;
    f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
end    
    
while abs((f/f_dash))>0.000001;
    z_f=z-f/f_dash;
    z=z_f;
    S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
    C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
    y=r1n+r2n+A*((z*S-1)/C^0.5);
    f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
    if z==0;
        f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
    else;
        f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
    end
end    

% Calculating the Lagrange functions and Velocities at r1 and r2

f=1-y/r1n;
g=(1/sqrt(u))*((y/C)^1.5*S)+A*sqrt(y)-(1/sqrt(u))*(y/C)^1.5;
f_dot=(sqrt(u)/r1n*r2n)*sqrt(y/C)*(z*S-1);
g_dot=1-(y/r2n);

v1_t=(r2-f*r1)/g;                                               % Transfer orbit velocity at r1
v2_t=(g_dot*r2-r1)/g;                                           % Transfer orbit velocity at r2
v1n_t=norm(v1_t);
v2n_t=norm(v2_t);

% Calculating velocities of Chaser at r1 and Target at r2
vpd_2=h_d/r2n;                                                  % Perpendicular component of velocity of the debri at 2
vrd_2=(u/h_d)*e*sind(ta_quarter);
vd_2=sqrt(vpd_2^2+vrd_2^2);
h_c=sqrt(u*r1n*(1+e));
vpc_1=h_c/r1n;                                                  % Perpedicular component of velocity of Chaser at 1 
vc_1=vpc_1;

% Calculating delta V required
delv1(1,1)=v1n_t-vc_1;
delv2(1,1)=vd_2-v2n_t;
delv(1,1)=abs(delv1(1,1))+abs(delv2(1,1));

%-------------------------------------------------------------------------------------------------------------------------------

% Now choosing a retrograde trajectory(i<90 deg). Rendezvous occurs at t=T/4 of debri orbit
if r1cr2(3,1)>=0;
    deltheta=360-acosd(dot(r1,r2)/r1r2);
else
    deltheta=acosd(dot(r1,r2)/r1r2);
end    

A=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));           

% Starting Newtons method to evaluate the z (z is  related to a and universal anomaly) 
z=0;                                                     % Initializing z
S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);  % S and C are Stumpff Functions 
C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
r1n=norm(r1);
r2n=norm(r2);
del_t=t_quarter;
y=r1n+r2n+A*((z*S-1)/C^0.5);
f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
if z==0;
    f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
else;
    f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
end    
    
while abs((f/f_dash))>0.01;
    z_f=z-f/f_dash;
    z=z_f;
    S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
    C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
    y=r1n+r2n+A*((z*S-1)/C^0.5);
    f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
    if z==0;
        f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
    else;
        f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
    end
end    

% Calculating the Lagrange functions and Velocities at r1 and r2

f=1-y/r1n;
g=(1/sqrt(u))*((y/C)^1.5*S)+A*sqrt(y)-(1/sqrt(u))*(y/C)^1.5;
f_dot=(sqrt(u)/r1n*r2n)*sqrt(y/C)*(z*S-1);
g_dot=1-(y/r2n);

v1_t=(r2-f*r1)/g;                                               % Transfer orbit velocity at r1
v2_t=(g_dot*r2-r1)/g;                                           % Transfer orbit velocity at r2
v1n_t=norm(v1_t);
v2n_t=norm(v2_t);

% Calculating delta V required
delv1(1,2)=v1n_t-vc_1;
delv2(1,2)=vd_2-v2n_t;
delv(1,2)=abs(delv1(1,2))+abs(delv2(1,2));

% Setting up Max DeltaV
if delv(1,2)<delv(1,1);
    delv_max=floor(delv(1,2));
else
    delv_max=floor(delv(1,1));
end    

%--------------------------------------------------------------------------------------------------------
% Checking Delta V for different True Anomaly in a Prograde Orbit

for i=1:80
    ta_chaser(1,i)=5*i;                                                    % Setting True Anomaly of chaser in degrees
    r(1,i)=(h_c^2)/(u*(1+e*cosd(ta_chaser(1,i))));
    vrc(1,i)=(u*e*sind(ta_chaser(1,i)))/h_c;                                 % Radial velocity of the chaser at that particular true anomaly   
    vpc(1,i)=h_c/r(1,i);                                                    % Perpendicular velocity of the chaser at that True Anomaly    
    vc(1,i)=sqrt(vpc(1,i)^2+vrc(1,i)^2);
    rc_cartesian=[r(1,i)*cosd(ta_chaser(1,i));r(1,i)*sind(ta_chaser(1,i));0];
    ia=95.2886;                                                             % Inclination Angle          
    DCM=[cosd(o),-sind(o),0;sind(o),cosd(o),0;0,0,1]*[1 0 0;0 cosd(ia) -sind(ia);0 sind(ia) cosd(ia)]*[cosd(w+ta_chaser(1,i)) -sind(w+ta_chaser(1,i)) 0; sind(w+ta_chaser(1,i)) cosd(w+ta_chaser(1,i)) 0; 0 0 1];
    r1(:,i)=DCM*rc_cartesian;                                    % Writing r1 in ECI frame

    % Choosing a prograde trajectory(i<90 deg). Rendezvous occurs at t=T/4 of debri orbit
    r1cr2=cross(r1(:,i),r2);
    r1r2=norm(r1(:,i))*norm(r2);
    if r1cr2(3,1)>=0;
        deltheta=acosd(dot(r1(:,i),r2)/r1r2);
    else
        deltheta=360-acosd(dot(r1(:,i),r2)/r1r2);
    end    
    A=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));  
    
    % Starting Newtons method to evaluate the z (z is  related to a and universal anomaly) 
    z=0;                                                     % Initializing z
    S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);  % S and C are Stumpff Functions 
    C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
    r1n=norm(r1(1,i));
    r2n=norm(r2);
    del_t=t_quarter;
    y=r1n+r2n+A*((z*S-1)/C^0.5);
    f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
    if z==0;
        f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
    else;
        f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
    end    
    while abs((f/f_dash))>0.01;                         % Executing Newtons Methods to find Z which is related to f and g  
        z_f=z-f/f_dash;
        z=z_f;
        S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
        C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
        y=r1n+r2n+A*((z*S-1)/C^0.5);
        f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
        if z==0;
            f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
        else;
            f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
        end
    end    
    % Calculating the Lagrange functions and Velocities at r1 and r2
    f=1-y/r1n;
    g=(1/sqrt(u))*((y/C)^1.5*S)+A*sqrt(y)-(1/sqrt(u))*(y/C)^1.5;
    f_dot=(sqrt(u)/r1n*r2n)*sqrt(y/C)*(z*S-1);
    g_dot=1-(y/r2n);
    v1_tp(:,i)=(r2-f*r1(:,i))/g;                                               % Transfer orbit velocity at r1
    v2_tp(:,i)=(g_dot*r2-r1(:,i))/g;                                           % Transfer orbit velocity at r2
    v1n_tp(1,i)=norm(v1_tp(:,i));
    v2n_tp(1,i)=norm(v2_tp(:,i));
    
    % Calculating delta V required
    delv1_pro(1,i)=v1n_tp(1,i)-vc(1,i);
    delv2_pro(1,i)=vd_2-v2n_tp(1,i);
    delv_pro(1,i)=abs(delv1_pro(1,i))+abs(delv2_pro(1,i));
end     

%--------------------------------------------------------------------------------------------------------
% Checking Delta V for different True Anomaly in a Retrograde Orbit

for i=1:80
    ta_chaser(1,i)=5*i;                                                    % Setting True Anomaly of chaser in degrees
    r(1,i)=(h_c^2)/(u*(1+e*cosd(ta_chaser(1,i))));
    vrc(1,i)=(u*e*sind(ta_chaser(1,i)))/h_c;                                 % Radial velocity of the chaser at that particular true anomaly   
    vpc(1,i)=h_c/r(1,i);                                                    % Perpendicular velocity of the chaser at that True Anomaly    
    vc(1,i)=sqrt(vpc(1,i)^2+vrc(1,i)^2);
    rc_cartesian=[r(1,i)*cosd(ta_chaser(1,i));r(1,i)*sind(ta_chaser(1,i));0];
    ia=95.2886;                                                             % Inclination Angle          
    DCM=[cosd(o),-sind(o),0;sind(o),cosd(o),0;0,0,1]*[1 0 0;0 cosd(ia) -sind(ia);0 sind(ia) cosd(ia)]*[cosd(w+ta_chaser(1,i)) -sind(w+ta_chaser(1,i)) 0; sind(w+ta_chaser(1,i)) cosd(w+ta_chaser(1,i)) 0; 0 0 1];
    r1(:,i)=DCM*rc_cartesian;                                    % Writing r1 in ECI frame

    % Choosing a Retrograde trajectory(i>90 deg). Rendezvous occurs at t=T/4 of debri orbit
    r1cr2=cross(r1(:,i),r2);
    r1r2=norm(r1(:,i))*norm(r2);
    if r1cr2(3,1)>=0;
        deltheta=360-acosd(dot(r1(:,i),r2)/r1r2);
    else
        deltheta=acosd(dot(r1(:,i),r2)/r1r2);
    end    
    A=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));  
    
    % Starting Newtons method to evaluate the z (z is  related to a and universal anomaly) 
    z=0;                                                     % Initializing z
    S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);  % S and C are Stumpff Functions 
    C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
    r1n=norm(r1(1,i));
    r2n=norm(r2);
    del_t=t_quarter;
    y=r1n+r2n+A*((z*S-1)/C^0.5);
    f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
    if z==0;
        f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
    else;
        f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
    end    
    while abs((f/f_dash))>0.01;                         % Executing Newtons Methods to find Z which is related to f and g  
        z_f=z-f/f_dash;
        z=z_f;
        S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
        C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
        y=r1n+r2n+A*((z*S-1)/C^0.5);
        f=((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*del_t;
        if z==0;
            f_dash=(sqrt(2)/40)*y^1.5+(A/8)*(sqrt(y)+A*sqrt(1/2*y));
        else;
            f_dash=((y/C)^1.5)*((1/2*z)*(C-1.5*S/C)+(3*S^2)/(4*C)) + (A/8)*((3*S*sqrt(y)/C)+A*sqrt(C/y));
        end
    end    
    % Calculating the Lagrange functions and Velocities at r1 and r2
    f=1-y/r1n;
    g=(1/sqrt(u))*((y/C)^1.5*S)+A*sqrt(y)-(1/sqrt(u))*(y/C)^1.5;
    f_dot=(sqrt(u)/r1n*r2n)*sqrt(y/C)*(z*S-1);
    g_dot=1-(y/r2n);
    v1_tr(:,i)=(r2-f*r1(:,i))/g;                                               % Transfer orbit velocity at r1
    v2_tr(:,i)=(g_dot*r2-r1(:,i))/g;                                           % Transfer orbit velocity at r2
    v1n_tr(1,i)=norm(v1_tr(:,i));
    v2n_tr(1,i)=norm(v2_tr(:,i));
    
    % Calculating delta V required
    delv1_retro(1,i)=v1n_tr(1,i)-vc(1,i);
    delv2_retro(1,i)=vd_2-v2n_tr(1,i);
    delv_retro(1,i)=abs(delv1_retro(1,i))+abs(delv2_retro(1,i));
end     

% Plotting Initial true anomaly Vs Delta V required for both trajectory

plot(ta_chaser,delv_retro,'r',ta_chaser,delv_pro,'b');
title('Initial True Anomaly of Chaser VS Delta V Required');
xlabel('Initial True Anomaly of Chaser (in degrees)');
ylabel('Net Delta V Required (in Km/sec)');
legend('Retrograde Orbit','Prograde Orbit');
