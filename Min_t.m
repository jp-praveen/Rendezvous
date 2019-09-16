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
a_chaser=a_debri+50;                                          % in km  
rp_debri=a_debri*(1-e);
rp_chaser=a_chaser*(1-e);
r1_cylindrical=rp_chaser

o=104.5636;                                             % RAAN
ia=95.2886;                                             % Inclination
w=189.5932;                                             % Argument of Periapsis

T_t=(2*pi*a_debri^(3/2))/u^0.5;
t_quarter=T_t/4;
T_c=(2*pi*a_chaser^(3/2))/u^0.5;
delv_max=18.8;

for i=0:52;
    twait=50*i;                                                   % Waiting time
    delV=50;
    ma_chaser=(2*pi*twait)/T_c;                     % ma=mean anomaly in radians
    if ma_chaser<pi;                                  % ea=eccentric anomaly     
        eai_chaser=ma_chaser+e/2;                    % initiating ea
    else
        eai_chaser=ma_chaser-e/2;
    end
    f=eai_chaser-e*sin(eai_chaser)-ma_chaser;
    f_dash=1-e*cos(eai_chaser);
    while abs((f/f_dash))>0.000001;                         % Newtons method to find eccentric anomaly 
        ea_chaser=eai_chaser-(f/f_dash);
        eai_chaser=ea_chaser;
        f=eai_chaser-e*sin(eai_chaser)-ma_chaser;
        f_dash=1-e*cos(eai_chaser); 
    end
    ta_chaser=2*atand(sqrt((1+e)/(1-e))*tand(ea_chaser/2))*180/3.14;    % ta=True Anomaly (in degrees) 
    r1_cylindrical=(h_d^2/u)*(1+e*cosd(ta_chaser))^(-1);
    r1_cartesian=[r1_cylindrical*cosd(ta_chaser);r1_cylindrical*sind(ta_chaser);0];
    DCMc=[cosd(o),-sind(o),0;sind(o),cosd(o),0;0,0,1]*[1 0 0;0 cosd(ia) -sind(ia);0 sind(ia) cosd(ia)]*[cosd(w+ta_chaser) -sind(w+ta_chaser) 0; sind(w+ta_chaser) cosd(w+ta_chaser) 0; 0 0 1];
    r1=DCMc*r1_cartesian;      
    for j=twait:50:t_quarter;
        j=dt;
        while delV>delV_max;
            t=twait+dt;
            ma_target=(2*pi*t)/T_t;                     % ma=mean anomaly in radians
            if ma_target<pi;                                  % ea=eccentric anomaly     
                eai_target=ma_target+e/2;                    % initiating ea
            else
                eai_target=ma_target-e/2;
            end
            f=eai_target-e*sin(eai_target)-ma_target;
            f_dash=1-e*cos(eai_target);
            while abs((f/f_dash))>0.000001;                         % Newtons method to find eccentric anomaly 
                ea_target=eai_target-(f/f_dash);
                eai_target=ea_target;
                f=eai_target-e*sin(eai_target)-ma_target;
                f_dash=1-e*cos(eai_target); 
            end
            ta_target=2*atand(sqrt((1+e)/(1-e))*tand(ea_target/2))*180/3.14;    % ta=True Anomaly (in degrees) 
            r2_cylindrical=(h_d^2/u)*(1+e*cosd(ta_target))^(-1);
            r2_cartesian=[r2_cylindrical*cosd(ta_target);r2_cylindrical*sind(ta_target);0];
            DCMt=[cosd(o),-sind(o),0;sind(o),cosd(o),0;0,0,1]*[1 0 0;0 cosd(ia) -sind(ia);0 sind(ia) cosd(ia)]*[cosd(w+ta_target) -sind(w+ta_target) 0; sind(w+ta_target) cosd(w+ta_target) 0; 0 0 1];
            r2=DCMt*r2_cartesian;
            
        
