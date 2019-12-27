% GIVEN R1,R2,dT FIND V1,V2 USING PRUSSINGS ALGORITHM
function [v1_short,v2_short,at_short,RAAN_short,inclination_short,perigee_short,v1_long,v2_long,at_long,RAAN_long,inclination_long,perigee_long] = lambert_prussing(r1,r2,transfer_time,dtheta,m)
u=132712440018;                             % in km^3*s^-2 

r1n=norm(r1);
r2n=norm(r2);
c=r2-r1;
cn=norm(c);
s1=(r1n+r2n+cn);
s2=(r1n+r2n-cn);
dt=transfer_time;

        

%SHORT PATH TRANSFER ORBIT
a0=100*r1n;
options = optimset('Display','off');
%F=sqrt((a_t(1)^3)/u)*((2*asin(sqrt(s1/(4*a_t)))-sin(2*asin(sqrt(s1/(4*a_t)))))-(2*asin(sqrt(s2/(4*a_t)))-sin(2*asin(sqrt(s2/(4*a_t))))))-dt;
%F_dash= -(a_t^3/u)^(1/2)*(s1/(4*a_t^2*(1 - s1/(4*a_t))^(1/2)*(s1/(4*a_t))^(1/2)) - s2/(4*a_t^2*(1 - s2/(4*a_t))^(1/2)*(s2/(4*a_t))^(1/2)) - (s1*cos(2*asin((s1/(4*a_t))^(1/2))))/(4*a_t^2*(1 - s1/(4*a_t))^(1/2)*(s1/(4*a_t))^(1/2)) + (s2*cos(2*asin((s2/(4*a_t))^(1/2))))/(4*a_t^2*(1 - s2/(4*a_t))^(1/2)*(s2/(4*a_t))^(1/2))) - (3*a_t^2*(sin(2*asin((s1/(4*a_t))^(1/2))) - sin(2*asin((s2/(4*a_t))^(1/2))) - 2*asin((s1/(4*a_t))^(1/2)) + 2*asin((s2/(4*a_t))^(1/2))))/(2*u*(a_t^3/u)^(1/2));
%F_dash=diff(F,a_t);
%dtheta=acos((r1n^2+r2n^2-cn^2)/(2*r1n*r2n));
sets1s2(s1,s2,u,dt,dtheta,m);
%at_short=fsolve(@solve_at_short,a0,options)                           % Semi major axis of the transfer orbit
if isnan(r1n)==0 & isnan(r2n)==0;
    at_short=solve_at_short_NR(a0);
else
    at_short=nan;
end
if at_short<0 & imag(at_short)==0;
    at_short=-1*at_short;
end
if imag(at_short)<0.000000001
    at_short=real(at_short);
end
if imag(at_short)==0 & abs(at_short)>s1/4 & abs(at_short)>s2/4 ;
    
    alpha_short=2*asin(sqrt(s1/(4*at_short)));
    beta_short=2*asin(sqrt(s2/(4*at_short)));

    g1_short=(1-r1n/at_short)*cot(alpha_short-beta_short)-(1-r2n/at_short)*csc(alpha_short-beta_short);
    g2_short=1-r1n/at_short;
       
    Eanomaly_trans_1_short=atan2(g1_short,g2_short);
    Eanomaly_trans_2_short=Eanomaly_trans_1_short+alpha_short-beta_short;

    e_short=(1-r1n/at_short)/cos(Eanomaly_trans_1_short);

    [RAAN_short,inclination_short,perigee_short,ta_1_short,ta_2_short]=orbitparameters(r1,r2,e_short,Eanomaly_trans_1_short,Eanomaly_trans_2_short);

    h_short= sqrt(r1n*u*(1+e_short*cosd(ta_1_short)));

    vp1_short=h_short/r1n;                                               % Transfer orbit velocity at r1
    vr1_short=(u/h_short)*e_short*sind(ta_1_short);
    vp2_short=h_short/r2n;                                               % Transfer orbit velocity at r2
    vr2_short=(u/h_short)*e_short*sind(ta_2_short);

    DCM_short=[cosd(RAAN_short),-sind(RAAN_short),0;sind(RAAN_short),cosd(RAAN_short),0;0,0,1]*[1 0 0;0 cosd(inclination_short) -sind(inclination_short);0 sind(inclination_short) cosd(inclination_short)]*[cosd(perigee_short+ta_1_short) -sind(perigee_short+ta_1_short) 0; sind(perigee_short+ta_1_short) cosd(perigee_short+ta_1_short) 0; 0 0 1];
    v1_short_1=[vr1_short;vp1_short;0];
    v1_short=DCM_short*v1_short_1;

    DCM_short=[cosd(RAAN_short),-sind(RAAN_short),0;sind(RAAN_short),cosd(RAAN_short),0;0,0,1]*[1 0 0;0 cosd(inclination_short) -sind(inclination_short);0 sind(inclination_short) cosd(inclination_short)]*[cosd(perigee_short+ta_2_short) -sind(perigee_short+ta_2_short) 0; sind(perigee_short+ta_2_short) cosd(perigee_short+ta_2_short) 0; 0 0 1];

    v2_short_1=[vr2_short;vp2_short;0];
    v2_short=DCM_short*v2_short_1;
else    
    v1_short=[nan;nan;nan];
    v2_short=[nan;nan;nan];
    RAAN_short=nan;
    inclination_short=nan;
    perigee_short=nan;
end

%LONG PATH TRANSFER ORBIT
a0=100*r1n;
options = optimset('Display','off');
%at_long=fsolve(@solve_at_long,a0,options);   % Semi major axis of the transfer orbit
if isnan(r1n)==0 & isnan(r2n)==0;
    at_long=solve_at_long_NR(a0);
else
    at_long=nan;
end
if at_long<0;
    at_long=-1*at_long;
%elseif imag(at_long)~=0
 %   at_long=real(at_long);
end
if imag(at_long)==0 & at_long>s1/4 & at_long>s2/4 & at_long>0;
    alpha_long=2*pi-2*asin(sqrt(s1/(4*at_long)));
    beta_long=2*asin(sqrt(s2/(4*at_long)));
 
    g1_long=(1-r1n/at_long)*cot(alpha_long-beta_long)-(1-r2n/at_long)*csc(alpha_long-beta_long);
    g2_long=1-r1n/at_long;
       
    Eanomaly_trans_1_long=atan2(g1_long,g2_long);
    Eanomaly_trans_2_long=Eanomaly_trans_1_long+alpha_long-beta_long;
 
    e_long=(1-r1n/at_long)/cos(Eanomaly_trans_1_long);
 
    [RAAN_long,inclination_long,perigee_long,ta_1_long,ta_2_long]=orbitparameters(r1,r2,e_long,Eanomaly_trans_1_long,Eanomaly_trans_2_long);
 
    h_long= sqrt(r1n*u*(1+e_long*cosd(ta_1_long)));
 
    vp1_long=h_long/r1n;                                               % Transfer orbit velocity at r1
    vr1_long=(u/h_long)*e_long*sind(ta_1_long);
    vp2_long=h_long/r2n;                                               % Transfer orbit velocity at r2
    vr2_long=(u/h_long)*e_long*sind(ta_2_long);
 
    DCM_long=[cosd(RAAN_long),-sind(RAAN_long),0;sind(RAAN_long),cosd(RAAN_long),0;0,0,1]*[1 0 0;0 cosd(inclination_long) -sind(inclination_long);0 sind(inclination_long) cosd(inclination_long)]*[cosd(perigee_long+ta_1_long) -sind(perigee_long+ta_1_long) 0; sind(perigee_long+ta_1_long) cosd(perigee_long+ta_1_long) 0; 0 0 1];
    v1_long_1=[vr1_long;vp1_long;0];
    v1_long=DCM_long*v1_long_1;
 
    DCM_long=[cosd(RAAN_long),-sind(RAAN_long),0;sind(RAAN_long),cosd(RAAN_long),0;0,0,1]*[1 0 0;0 cosd(inclination_long) -sind(inclination_long);0 sind(inclination_long) cosd(inclination_long)]*[cosd(perigee_long+ta_2_long) -sind(perigee_long+ta_2_long) 0; sind(perigee_long+ta_2_long) cosd(perigee_long+ta_2_long) 0; 0 0 1];
 
    v2_long_1=[vr2_long;vp2_long;0];
    v2_long=DCM_long*v2_long_1;
else
    v1_long=[nan;nan;nan];
    v2_long=[nan;nan;nan];
    RAAN_long=nan;
    inclination_long=nan;
    perigee_long=nan;
end


