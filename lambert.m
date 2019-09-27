% GIVEN INITAL POSITIONS AND TRANSFER TIME FIND THE VELOCITIES OF THE PROGRADE TRANSFER ORBIT
function [v1,v2] = lambert(r1,r2,transfer_time)
u=398588.738;                             % in km^3*s^-2 

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
