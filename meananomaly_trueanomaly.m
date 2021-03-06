% GIVEN MEAN ANOMALY FIND THE EQUIVALENT TRUE ANOMALY IN THE ORBIT
function [true_anomaly] = meananomaly_trueanomaly(mean_anomaly_deg)
mat=getparameters;
e=mat(1,6);
mean_anomaly=mean_anomaly_deg*3.14/180;
setmeananomaly(mean_anomaly);
if mean_anomaly<pi;                                             % initializing eccentric anomaly     
    eccentric_anomaly1=mean_anomaly+e/2;                                % initiating ea
else
    eccentric_anomaly1=mean_anomaly-e/2;
end
eccentric_anomaly=fzero(@solve_eccentricanomaly,eccentric_anomaly1);

true_anomaly_rad=2*atan(sqrt((1+e)/(1-e))*tan(eccentric_anomaly/2));    
true_anomaly=true_anomaly_rad*180/3.14;                      % ta=True Anomaly (in degrees) 
