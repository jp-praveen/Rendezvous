% If Approach==1 ---> Greedy only
% If Approach==other than 1 ---> Greedy Heuristic
function [sequence_best,dv,time_elapsed] = greedy(totalnode,startnode,dvmax,approach) 
u=398588.738;                             % in km^3*s^-2 
tn=totalnode;
s=startnode;
dv=0;
sequence_best=[s];
time_elapsed=0;
twait=0;

fname = 'data_check_paper_1to50.txt';

% Open the TLE file and read TLE elements
fid = fopen(fname, 'r');

for i=1:2*tn;
    % read first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    Cnum(1,i) = str2num(tline(3:7));      			        % Catalog Number (NORAD)
    epoch(1,i) = str2num(tline(19:32));              % Epoch
    Etype(1,i) = str2num(tline(63));                          % Ephemeris Type
    Enum(1,i)  = str2num(tline(65:end));             % Element Number
    
    % read second line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    inclination(1,i) = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    RAAN(1,i) = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e(1,i) = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    perigee(1,i) = str2num(tline(35:42));              % Argument of Perigee (degrees)
    mean_anomaly(1,i) = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    no(1,i) = str2num(tline(53:63));                 % Mean Motion
    a(1,i) =( u/(no(1,i)*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    rNo(1,i) = str2num(tline(64:68));                % Revolution Number at Epoch
end
fclose(fid);

[x,y,z,xdot,ydot,zdot]=test_sgp4;

for i=1:tn;
        debri_position_1(:,i)=[x(i,1);y(i,1);z(i,1)];          % Initial position of debris; nth column implies nth debri data
        debri_velocity_1(:,i)=[xdot(i,1);ydot(i,1);zdot(i,1)]; % Initial velocity of debris; nth column implies nth debri data
        h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));                      % specific angluar momentum   
        T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);                  % Time Period
        t_initial(1,i)=(mean_anomaly(1,i)*T(1,i))/(2*180);
        true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));     
end
%setdebrisdata(debri_position_1,debri_velocity_1,h,T,t_initial,true_anomaly);

while dv<dvmax;
    [time_elapsed,dv_current,next_node] = realistic_position(tn,s,twait,sequence_best,dvmax,approach,debri_position_1,debri_velocity_1,h,T,t_initial,true_anomaly,e);
    if dv+dv_current<dvmax;
        s=next_node;
        twait=twait+time_elapsed;
        dv=dv+dv_current
        sequence_best=[sequence_best s];
    else
        break;
    end
end
time_elapsed=twait;
  