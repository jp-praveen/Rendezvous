% GIVEN TRUE ANOMALY IN AN ORBIT FIND THE EQUIVALENT POSITION VECTOR IN ECI FRAME
function [r] = tanomaly_rad(true_anomaly,h)
mat=getparameters;
u=mat(1,1);
e=mat(1,6);
RAAN=mat(1,2);
inclination=mat(1,3);
perigee=mat(1,4);
r_cylindrical=(h^2/u)*(1+e*cosd(true_anomaly))^(-1);
r_cartesian=[r_cylindrical*cosd(true_anomaly);r_cylindrical*sind(true_anomaly);0];
DCM=[cosd(RAAN),-sind(RAAN),0;sind(RAAN),cosd(RAAN),0;0,0,1]*[1 0 0;0 cosd(inclination) -sind(inclination);0 sind(inclination) cosd(inclination)]*[cosd(perigee) -sind(perigee) 0; sind(perigee) cosd(perigee) 0; 0 0 1];
r=DCM*r_cartesian;
