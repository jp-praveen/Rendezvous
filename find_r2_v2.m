% Given r1 and v1 in an orbit, find r2 and v2 at a different time
% [Algorithm 3.4 in Orbital Mechanics by Howard Curtis
function [r2,v2,alpha,universal_anomaly] = find_r2_v2(r1,v1,dt,T)
u=398588.738;                             % in km^3*s^-2 
r0=sqrt(dot(r1,r1));
v0=sqrt(dot(v1,v1));
v_r=dot(r1,v1)/r0;                 % Radial velocity
alpha=(2/r0)-(v0^2/u);
if dt>T;
    ratio=floor(dt/T);
    dt=dt-ratio*T;
end

universal_anomaly_estimate=sqrt(u)*abs(alpha)*dt;
z=alpha*universal_anomaly_estimate^2;
S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
f_dash=(r0*v_r*universal_anomaly_estimate/sqrt(u))*(1-alpha*universal_anomaly_estimate^2*S)+(1-alpha*r0)*universal_anomaly_estimate^2*C+r0;   
f=((r0*v_r/sqrt(u))*universal_anomaly_estimate^2*C)+((1-alpha*r0)*S*universal_anomaly_estimate^3)+r0*universal_anomaly_estimate-sqrt(u)*dt;
%fun=f/f_dash;
%universal_anomaly=fzero(@solve_r2_v2,universal_anomaly_estimate1,alpha,r0,v0,v_r,u,dt);

universal_anomaly=universal_anomaly_estimate-f/f_dash;
p=0;
while abs(f/f_dash)>0.000001;
    p=p+1;
    universal_anomaly=universal_anomaly_estimate-f/f_dash;
    universal_anomaly_estimate=universal_anomaly;
    z=alpha*universal_anomaly_estimate^2;
    S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
    C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
    f_dash=(r0*v_r*universal_anomaly_estimate/sqrt(u))*(1-alpha*universal_anomaly_estimate^2*S)+(1-alpha*r0)*universal_anomaly_estimate^2*C+r0;   
    f=((r0*v_r/sqrt(u))*universal_anomaly_estimate^2*C)+((1-alpha*r0)*S*universal_anomaly_estimate^3)+r0*universal_anomaly_estimate-sqrt(u)*dt;
    
    if p>10000;
        break;
    end
end
check=f/f_dash;
if p<10000;
    f=1-universal_anomaly^2*C/r0;
    g=dt-universal_anomaly^3*S/sqrt(u);

    r2=f*r1+g*v1;
    r2n=norm(r2);

    f_dot=(sqrt(u)/(r2n*r0))*(alpha*S*universal_anomaly^3-universal_anomaly);
    g_dot=1-universal_anomaly^2*C/r2n;

    v2=f_dot*r1+g_dot*v1;
else
    r2=nan;
    v2=nan;
end

    