function F = solve_r2_v2(universal_anomaly_estimate,alpha,r0,v0,v_r,u,dt)
z=alpha*universal_anomaly_estimate^2;
S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
f_dash=(r0*v_r*universal_anomaly_estimate/sqrt(u))*(1-alpha*universal_anomaly_estimate^2*S)+(1-alpha*r0)*universal_anomaly_estimate^2*C+r0;   
f=((r0*v_r/sqrt(u))*universal_anomaly_estimate^2*C)+((1-alpha*r0)*S*universal_anomaly_estimate^3)+r0*universal_anomaly_estimate-sqrt(u)*dt;
F(1)=f/f_dash;