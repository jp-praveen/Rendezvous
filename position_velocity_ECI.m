% GIVEN POSITION AND VELOCITY VECTOR AT A GIVEN TRUE ANOMALY, FIND POSITION VECTOR IN ECI FRAME
function [r,v] = position_velocity_ECI (position_vector,velocity_vector,true_anomaly)
mat=getparameters;
RAAN=mat(1,2);
inclination=mat(1,3);
perigee=mat(1,4);
DCM=[cosd(RAAN),-sind(RAAN),0;sind(RAAN),cosd(RAAN),0;0,0,1]*[1 0 0;0 cosd(inclination) -sind(inclination);0 sind(inclination) cosd(inclination)]*[cosd(perigee+true_anomaly) -sind(perigee+true_anomaly) 0; sind(perigee+true_anomaly) cosd(perigee+true_anomaly) 0; 0 0 1];
r=DCM*transpose(position_vector);
v=DCM*transpose(velocity_vector);
