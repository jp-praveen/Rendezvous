% GIVEN INITAL POSITIONS AND TRANSFER TIME FIND THE VELOCITIES OF THE PROGRADE TRANSFER ORBIT
function [v1,v2,RAAN,inclination,perigee,true_anomaly_1,true_anomaly_2] = lambert(r1,r2,transfer_time)
u=398588.738;                             % in km^3*s^ 

r1cr2=cross(r1,r2);
r1r2=norm(r1)*norm(r2);
if r1cr2(3,1)>=0;
    deltheta=acosd(dot(r1,r2)/r1r2);                                  % deltheta is change in true anomaly in transfer orbit or Transfer angle
else
    deltheta=360-acosd(dot(r1,r2)/r1r2);                              % deltheta is in degrees  
end    
r1n=norm(r1);
r2n=norm(r2);
dt=transfer_time;
A=sind(deltheta)*sqrt(r1r2/(1-cosd(deltheta)));
setsolveparameters(u,dt,r1n,r2n,A);
z0=[0];
z=fzero(@solve_torbit,z0);

S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);  % S and C are Stumpff Functions 
C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
y=r1n+r2n+A*((z*S-1)/C^0.5);

% Calculating the Lagrange functions and Velocities at r1 and r2

f=1-y/r1n;
g=(A*sqrt(y)/sqrt(u));
f_dot=(sqrt(u)/r1n*r2n)*sqrt(y/C)*(z*S-1);
g_dot=1-(y/r2n);
v1=(r2-f*r1)/g;                                               % Transfer orbit velocity at r1
v2=(g_dot*r2-r1)/g;                                           % Transfer orbit velocity at r2

% Calculating the Orbital Elements of the transfer orbit
evec = ((norm(v1)^2-u/norm(r1))*r1-dot(r1,v1)*v1)/u;
e = norm(evec);
v1_r=dot(r1,v1)/r1n;
if v1_r>=0
    true_anomaly_1=acosd(dot(evec,r1)/(e*r1n));
else
    true_anomaly_1=360-acosd(dot(evec,r1)/(e*r1n));
end
v2_r=dot(r2,v2)/r2n;
if v2_r>=0
    true_anomaly_2=acosd(dot(evec,r2)/(e*r2n));
else
    true_anomaly_2=360-acosd(dot(evec,r2)/(e*r2n));
end
h=cross(transpose(r1),transpose(v1));
K=[0,0,1];
n=cross(K,h);
inclination=acosd((h(1,3)/norm(h)));
RAAN = acosd(n(1,1)/norm(n));
if n(1,2)<0;
   RAAN = 360-RAAN;
end

perigee = acosd(dot(n,evec)/(norm(n)*e));

if evec(3,1)<0;
    perigee = 360-perigee;
end

