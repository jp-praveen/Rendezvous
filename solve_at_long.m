function F = solve_at_long(a_t)
matrix=gets1s2;
s1=matrix(1,1);
s2=matrix(1,2);
u=matrix(1,3);
dt=matrix(1,4);
alpha=2*pi-2*asin(sqrt(s1/(4*a_t)));                   
beta=2*asin(sqrt(s2/(4*a_t)));

F(1)=sqrt((a_t(1)^3)/u)*((alpha-sin(alpha))-(beta-sin(beta)))-dt;
