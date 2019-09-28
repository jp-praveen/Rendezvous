% GIVEN R1,R2,dT FIND V1,V2 USING PRUSSINGS ALGORITHM
function [v1,v2,RAAN_trans,inclination_trans,perigee_trans] = lambert(r1,r2,transfer_time)
u=398588.738;                             % in km^3*s^-2 

r1n=norm(r1);
r2n=norm(r2);
c=r2-r1;
cn=norm(c);
s1=(r1n+r2n+cn);
s2=(r1n+r2n-cn);
dt=transfer_time;
sets1s2(s1,s2,u,dt);
        
a0=100*r1n;
at1=fzero(@solve_at,a0);                            % Semi major axis of the transfer orbit
at2=fsolve(@solve_at,a0);
at=(at1+at2)*0.5;
alpha=2*asin(sqrt(s1/(4*at)));
beta=2*asin(sqrt(s2/(4*at)));

g1=(1-r1n/at)*cot(alpha-beta)-(1-r2n/at)*csc(alpha-beta);
g2=1-r1n/at;
       
Eanomaly_trans_1=atan2(g1,g2);
Eanomaly_trans_2=Eanomaly_trans_1+alpha-beta;
e_t=(1-r1n/at)/cos(Eanomaly_trans_1);
[RAAN_trans,inclination_trans,perigee_trans,ta_1,ta_2]=orbitparameters(r1,r2,e_t,Eanomaly_trans_1,Eanomaly_trans_2);

h_trans= sqrt(r1n*u*(1+e_t*cosd(ta_1)));
vp1_transfer=h_trans/r1n;                                               % Transfer orbit velocity at r1
vr1_transfer=(u/h_trans)*e_t*sind(ta_1);
vp2_transfer=h_trans/r2n;                                               % Transfer orbit velocity at r2
vr2_transfer=(u/h_trans)*e_t*sind(ta_2);

DCM_trans=[cosd(RAAN_trans),-sind(RAAN_trans),0;sind(RAAN_trans),cosd(RAAN_trans),0;0,0,1]*[1 0 0;0 cosd(inclination_trans) -sind(inclination_trans);0 sind(inclination_trans) cosd(inclination_trans)]*[cosd(perigee_trans+ta_1) -sind(perigee_trans+ta_1) 0; sind(perigee_trans+ta_1) cosd(perigee_trans+ta_1) 0; 0 0 1];
v1_transfer=[vr1_transfer;vp1_transfer;0];
v1=DCM_trans*v1_transfer;

DCM_trans=[cosd(RAAN_trans),-sind(RAAN_trans),0;sind(RAAN_trans),cosd(RAAN_trans),0;0,0,1]*[1 0 0;0 cosd(inclination_trans) -sind(inclination_trans);0 sind(inclination_trans) cosd(inclination_trans)]*[cosd(perigee_trans+ta_2) -sind(perigee_trans+ta_2) 0; sind(perigee_trans+ta_2) cosd(perigee_trans+ta_2) 0; 0 0 1];

v2_transfer=[vr2_transfer;vp2_transfer;0];
v2=DCM_trans*v2_transfer;
