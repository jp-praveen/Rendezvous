% Solving for eccentric anomaly using Fsolve
function F = solve_eccentricanomaly(eccentric_anomaly)
mat=getparameters;
e=mat(1,6);
matrix=getmeananomaly;
mean_anomaly=matrix;
F(1)=eccentric_anomaly-e*sin(eccentric_anomaly)-mean_anomaly;
