% If Approach==1 ---> Greedy only
% If Approach==2 ---> Greedy Heuristic with heuristic_cost(i)=dvmax/depth search sequence_length;
% If Approach==3 ---> Greedy Heuristic with heuristic_cost(i)=dvmax/outdegree(i);
% If Approach==4 ---> Greedy Heuristic with heuristic_cost(i)=dvmax/depth search sequence_length with dvmax_available;


function [sequence_best,dv,time_elapsed] = greedy(totalnode,startnode,dvmax,approach,previous_sequence) 
u=132712440018;                             % in km^3*s^-2 
tn=totalnode;
s=startnode;
dv=0;
sequence_best=[s];
time_elapsed=0;
twait=0;

%fname = 'data_check_paper_1to50.txt';

% Open the TLE file and read TLE elements
%fid = fopen(fname, 'r');

%for i=1:2*tn;
    % read first line
 %   tline = fgetl(fid);
  %  if ~ischar(tline)
   %     break
    %end
    %Cnum(1,i) = str2num(tline(3:7));      			        % Catalog Number (NORAD)
    %epoch(1,i) = str2num(tline(19:32));              % Epoch
    %Etype(1,i) = str2num(tline(63));                          % Ephemeris Type
    %Enum(1,i)  = str2num(tline(65:end));             % Element Number
    
    % read second line
    %tline = fgetl(fid);
    %if ~ischar(tline)
     %   break
    %end
    %inclination(1,i) = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    %RAAN(1,i) = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    %e(1,i) = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    %perigee(1,i) = str2num(tline(35:42));              % Argument of Perigee (degrees)
    %mean_anomaly(1,i) = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    %no(1,i) = str2num(tline(53:63));                 % Mean Motion
    %a(1,i) =( u/(no(1,i)*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    %rNo(1,i) = str2num(tline(64:68));                % Revolution Number at Epoch
%end
%fclose(fid);
%------------------------------------------------------------------------------

% CHECK FOR SAME FILE NAME 
%[x,y,z,xdot,ydot,zdot]=test_sgp4;

%------------------------------------------------------------------------------
%for i=1:tn;
 %       debri_position_1(:,i)=[x(i,1);y(i,1);z(i,1)];          % Initial position of debris; nth column implies nth debri data
  %      debri_velocity_1(:,i)=[xdot(i,1);ydot(i,1);zdot(i,1)]; % Initial velocity of debris; nth column implies nth debri data
   %     h(1,i)=sqrt(a(1,i)*u*(1-e(1,i)^2));                      % specific angluar momentum   
    %    T(1,i)=((u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5))^(-1)*(2*pi);                  % Time Period
     %   t_initial(1,i)=(mean_anomaly(1,i)*T(1,i))/(2*180);
      %  true_anomaly(1,i)=ma_ta(mean_anomaly(1,i),e(1,i));     
%end
%setdebrisdata(debri_position_1,debri_velocity_1,h,T,t_initial,true_anomaly);
au=1.496e+8;
debri_position_1=[-0.17206*au,-1.02848*au,-0.34056*au,0.28759*au,-0.90624*au; 0.96812*au,-0.26939*au,-1.00392*au,1.00171*au,-0.54981*au; -0.00001*au,0.02359*au,0.00388*au,0.00877*au,-0.01290*au];
debri_velocity_1=[-29.815,9.136,27.819,-26.892,16.472;-5.325,-27.02,-7.605,10.019,-21.856;0.0,0.029,-0.060,-0.659,-0.056];
%debri_position_1=[-0.17206,0.94499;0.96812,-0.21205;-0.00001,0.00322]*au;
%debri_velocity_1=[-29.815,7.668;-5.325,30.242;0,0.102];
for i=1:tn
    radius_vec=debri_position_1(:,i);
    velocity_vec=debri_velocity_1(:,i)
    norm_r=norm(radius_vec);
    norm_v=norm(velocity_vec);

    e(:,i)=((norm_v^2-(u/norm_r))*radius_vec-(dot(radius_vec,velocity_vec)*velocity_vec))/u;
    e_norm(i)=norm(e(:,i));
    E_chaser=norm_v^2/2 - u/norm_r;
    a(i)=-u/(2*E_chaser);
    a_au=a(i)/au;
    if dot(radius_vec,velocity_vec)>0
        true_anomaly(i)=acos(dot(e(:,i),radius_vec)/(e_norm(i)*norm_r));
        true_anomaly(i)=true_anomaly(i)*180/pi;
    else
        true_anomaly(i)=2*pi-acos(dot(e(:,i),radius_vec)/(e_norm(i)*norm_r));
        true_anomaly(i)=true_anomaly(i)*180/pi;
    end
    ecc_anomaly(i)=2*atan(sqrt((1-e_norm(i))/(1+e_norm(i)))*tand(true_anomaly(i)/2));
    mean_anomaly(i)=ecc_anomaly(i)-e_norm(i)*sin(ecc_anomaly(i));
    mean_anomaly(i)=mean_anomaly(i)*180/pi;
    if mean_anomaly(i)<0
        mean_anomaly(i)=360+mean_anomaly(i);
    end
end


%inclination=[1.2,1.2798,0.23523,1.43929,0.72205,1.26705,0.55071,0.26337,1.40948,3.57355,4.80881,3.10322,1.69487,2.62074,0.77337,6.32151,2.20908];
%e=[e_chaser_norm,0.074019,0.0604714,0.083094,0.1448353,0.0515307,0.1065636,0.0684246,0.0985837,0.2128808,0.2194123,0.0862374,0.0967254,0.1109538,0.0838748,0.1698532,0.1389606];
%a=[a_chaser_au,1.0377497,1.0537339,1.0099643,0.9593424,1.0617043,1.0382537,1.033219,1.0235198,1.0125947,1.0324714,1.0621126,1.1441219,0.9115237,0.9936781,1.0323256,1.0075469];
%a=a*au;

%mean_anomaly=[mean_anomaly,220.55456,242.72351,57.934,146.08995];
%for i=1:5
 %  if true_anomaly(i)<0
  %      true_anomaly(i)=360+true_anomaly(i);
   % end
%end
mean_anomaly=mean_anomaly*pi/180;
for i=1:tn;
    h(1,i)=sqrt(a(1,i)*u*(1-e_norm(1,i)^2));                      % specific angluar momentum   
    T(1,i)=((u^2/h(1,i)^3)*(1-e_norm(1,i)^2)^(1.5))^(-1)*(2*pi);
    t_initial(1,i)=mean_anomaly(1,i)*T(1,i)/(2*pi);
end
mean_anomaly=mean_anomaly*180/pi;   
for twait=2026.5;
    twait
    twait=twait*86400;
    dv=0;
    dv_current=0;
    s=startnode;
    sequence_best=[s];
    
    while dv<dvmax & length(sequence_best)<tn
        %twait
        dvmax_available=dvmax-dv;
        [time_elapsed,dv_current,next_node] = realistic_position(tn,s,twait,sequence_best,previous_sequence,dvmax,dvmax_available,approach,debri_position_1,debri_velocity_1,h,T,t_initial,true_anomaly,e);
        %time_elapsed
        if dv+dv_current<dvmax;
            s=next_node;
            twait=twait+time_elapsed;
            dv=dv+dv_current
            sequence_best=[sequence_best s]
        else
            break;
        end
    end
end
time_elapsed=twait;
  