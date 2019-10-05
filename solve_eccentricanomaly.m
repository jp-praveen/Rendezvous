function F = solve_eccentricanomaly(eccentric_anomaly)
%mat=getparameters;
%e=mat(1,6);
matrix(:,1)=getmeananomaly;
mean_anomaly=matrix(1,1);
e=matrix(2,1);

F(1)=eccentric_anomaly-e*sin(eccentric_anomaly)-mean_anomaly;
