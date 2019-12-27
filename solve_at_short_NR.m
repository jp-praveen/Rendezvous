function [a1] = solve_at_short_NR(a0);
syms alp(a_t) bet(a_t) F(a_t) F_dash(a_t);
matrix=gets1s2;

s1=matrix(1,1);
s2=matrix(1,2);
u=matrix(1,3);
dt=matrix(1,4);
dtheta=matrix(1,5);
m=matrix(1,6);
if dtheta>1800;
    alp(a_t)=2*pi-2*asin(sqrt(s1/(4*a_t)));                   
    bet(a_t)=-2*asin(sqrt(s2/(4*a_t)));
    F(a_t)=sqrt((a_t(1)^3)/u)*(2*m*pi+((2*pi-2*asin(sqrt(s1/(4*a_t))))-sin(2*pi-2*asin(sqrt(s1/(4*a_t)))))-((-2*asin(sqrt(s2/(4*a_t))))-sin(-2*asin(sqrt(s2/(4*a_t))))))-dt;
    F_dash(a_t)=diff(F); 
    F_ddash(a_t)=diff(F_dash);
    F_ddash(a_t)=diff(F_dash);
else
    alp(a_t)=2*asin(sqrt(s1/(4*a_t)));                   
    bet(a_t)=2*asin(sqrt(s2/(4*a_t)));
    F(a_t)=sqrt((a_t(1)^3)/u)*(2*m*pi+((2*asin(sqrt(s1/(4*a_t))))-sin(2*asin(sqrt(s1/(4*a_t)))))-((2*asin(sqrt(s2/(4*a_t))))-sin(2*asin(sqrt(s2/(4*a_t))))))-dt;
    F_dash(a_t)=diff(F); 
    F_ddash(a_t)=diff(F_dash);
end


i=0;
%a_t=a0;
%subs(F)
%subs(F_dash)
%eval(F(a0))
%eval(F_dash(a0))

while abs((eval(F(a0))/eval(F_dash(a0))))>0.0000001
    abs((eval(F(a0))/eval(F_dash(a0))));
    i=i+1;
    %if eval(F_dash(a0))~=0;
     %   a1=1-eval(F(a0))/eval(F_dash(a0));
    %end
    a1=a0-(2*eval(F(a0))*eval(F_dash(a0)))/((2*eval(F_dash(a0))^2)-eval(F(a0))*eval(F_ddash(a0)));
    a0=a1;
    %a_t=a0;
    if i>50;
        break;
    end
end
%a1
%abs((eval(F(a0))/eval(F_dash(a0))))
    