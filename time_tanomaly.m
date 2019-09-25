% GIVEN TIME 't' IN ORBIT FIND TRUE ANOMALY
function [true_anomaly] = time_tanomaly(t,a,e)
u=398588.738;                                                % in km^3*s^-2  
T=(2*pi*a^(3/2))/sqrt(u);
mean_anomaly=(2*pi*t)/T;                                     % mean anomaly in radians 
mean_anomaly_deg=mean_anomaly*180/3.14;                      % mean anomaly in degrees 
true_anomaly=ma_ta(mean_anomaly_deg);                        % ta=True Anomaly (in degrees) 
