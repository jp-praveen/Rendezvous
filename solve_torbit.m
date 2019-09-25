% Solving for z using Fsolve function
function F = solve_torbit(z)
matrix=getsolveparameters; 
u=matrix(1,1);
dt=matrix(1,2);
r1n=matrix(1,3);
r2n=matrix(1,4);
A=matrix(1,5);
S=(1/6)-(z/120)+(z^2/5040)-(z^3/362880)+(z^4/39916800);   
C=(1/2)-(z/24)+(z^2/720)-(z^3/40320)+(z^4/3628800);
y=r1n+r2n+A*((z*S-1)/C^0.5);
if dt~=0;
    F(1) = ((y/C)^1.5)*S+A*(y^0.5)-(u^0.5)*dt;
else
    F(1) = ((y/C)^1.5)*S+A*(y^0.5);
end
