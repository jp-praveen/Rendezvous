function [RAAN,inclination,perigee,Tanomaly_1] = orbitparameters(r1,r2,e,Eanomaly_1,Eanomaly_2)
Tanomaly_1=2*atand(sqrt((1+e)/(1-e))*tand(Eanomaly_1*0.5))*180/3.14;
if Tanomaly_1<0;
    Tanomaly_1=360+Tanomaly_1;
end
Tanomaly_2=2*atand(sqrt((1+e)/(1-e))*tand(Eanomaly_2*0.5))*180/3.14
del_theta=Tanomaly_2-Tanomaly_1;
if del_theta<180;
    H=cross(r1,r2);
else
    H=cross(r2,r1);
end
inclination=atand(sqrt(H(1,1)^2+H(2,1)^2)/H(3,1));
RAAN=atand(H(1,1)/-H(2,1));

u_t1=atan2d((-r1(1,1)*cosd(inclination)*sind(RAAN)+r1(2,1)*cosd(inclination)*cosd(RAAN)+r1(3,1)*sind(inclination)),(r1(1,1)*cosd(RAAN)+r1(2,1)*sind(RAAN)));
perigee=(u_t1-Tanomaly_1);
if perigee<0;
    perigee=360+perigee;
end
