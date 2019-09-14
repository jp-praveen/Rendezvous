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
            542.13411088    -3934.00754565    -5496.24138485 
             XDOT             YDOT             ZDOT    [km/s]
            -2.04410472        5.83405751       -4.77434023

e=0.572671;
u=398588.738                                                 %in km^3*s^-2                
pos_d=[542.13411088 -3934.00754565 -5496.24138485];          %in km
vel_d=[-2.04410472 5.83405751 -4.77434023];                  %in kmpersec
r_d=norm(pos_d);
v_d=norm(vel_d);
vr_d=dot(vel_d,pos_d)/r_d;
vp_d=sqrt(v_d^2-vr_d^2);
h_d=r_d*vp_d;
a_debri=(h_d^2/u)/(1-e^2);                                     %in km
a_target=a_debri+50                                          %in km  
rp_debri=a_debri*(1-e);
rp_target=a_target*(1-e);

%Assuming both target and chaser are at periapsis at t=0

r1=rp_target;
r2=rp_debri;
