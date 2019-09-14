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

ta_quarter=2*atan(sqrt((1+e)/(1-e))*tan(ea_quarter/2))*180/3.14;    % ta=True Anomaly (in degrees) 
r_quarter=(h_d^2/u)*(1+e*cos(ta_quarter))^(-1);
r2_cylindrical=r_quarter;

% Converting radial vectors to ECI frame
ta_debri=ta_quarter;
ta_chaser=0;                                             % True Anomaly of chaser
r1_cartesian=[r1_cylindrical*cos(ta_chaser);r1_cylindrical*sin(ta_chaser);0];
r2_cartesian=[r2_cylindrical*cos(ta_debri);r1_cylindrical*sin(ta_debri);0];

o=104.5636;                                             % RAAN
i=95.2886;                                              % Inclination
w=189.5932;                                             % Argument of Periapsis
DCM=[cos(o),-sin(o),0;sin(o),cos(o),0;0,0,1]*[1 0 0;0 cos(i) -sin(i);0 sin(i) cos(i)]*[cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];

r1=DCM*r1_cartesian;                                    % Writing r1 and r2 in ECI frame
r2=DCM*r2_cartesian;

%Assuming both target and chaser are at periapsis at t=0 and choosing a prograde trajectory(i<90 deg). Rendezvous occurs at t=T/4 of %debri orbit
r1cr2=cross(r1,r2);
r1r2=norm(r1)*norm(r2);
if r1cr2(3,1)>=0;
    deltheta=acos(dot(r1,r2)/r1r2);
else
    deltheta=360-acos(dot(r1,r2)/r1r2);

A=sin(deltheta)*sqrt(r1r2/(1-cos(deltheta)));





