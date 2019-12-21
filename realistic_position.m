function [time_elapsed,dv,next_node] = realistic_position(totalnode,startnode,twait,sequence_best,previous_sequence,dvmax,dvmax_available,approach,debri_position_1,debri_velocity_1,h,T,t_initial,true_anomaly,e)        
tn=totalnode;
s=startnode;
method=3;
u=398588.738;                             % in km^3*s^-2 
                                
for i=1:tn
    r=debri_position_1(:,i);
    v=debri_velocity_1(:,i);
    [r1,v1,alpha,universal_anomaly]=find_r2_v2(r,v,twait,T(1,i));
    debri_position_new(i,:)=r1;
    debri_velocity_new(i,:)=v1;
end

r1=debri_position_new(s,:);
v1=debri_velocity_new(s,:);
r1=transpose(r1);
v1=transpose(v1);
ma_deg_chaser=(u^2/h(s)^3)*(1-e(s)^2)^(1.5)*(t_initial(s)+twait)*180/pi;
t_initial(s);
%twait
%h(s)
%e(s)
ta_chaser=ma_ta(ma_deg_chaser,e(s));

for i=1:tn;
    
    if ismember(i,sequence_best)==0 & ismember(i,previous_sequence)==0;
        r2_ini=debri_position_new(i,:);
        v2_ini=debri_velocity_new(i,:);
        k=0;
        if approach==2;
            sequence_heuristic = depth_search(50,i,dvmax);
            sequence_length=length(sequence_heuristic);
            if sequence_length~=0;
                heuristic_cost(i)=dvmax/sequence_length;
            else
                heuristic_cost(i)=1000;
            
            end
        elseif approach==3
            %[possible_start,possible_end]=depth_search.sf1_depth_search(totalnode);
            %[possible_start]=depth_search('sf1_depth_search',totalnode);
            load matrixmat.mat;
            matrixmat=matrixmat{:,:};
            matrixmat(1,:) = [];                    % Since 1st row and column contains the serial numbers
            matrixmat(:,1) = [];
            possible_start=[];
            possible_end=[];
            tn=totalnode;
            for j1=1:tn;
                for k1=1:50;
                    if isnan(matrixmat(j1,k1))==0;
                        possible_start=[possible_start j1];
                        possible_end=[possible_end k1];
                    end
                end
            end
            connections=possible_start(possible_start==s);
            outdegree=length(connections);
            if outdegree~=0;
                heuristic_cost(i)=dvmax/outdegree;
            else
                heuristic_cost(i)=1000;
            end
        elseif approach==4;
            sequence_heuristic = depth_search(50,i,dvmax_available);
            sequence_length=length(sequence_heuristic);
            if sequence_length~=0;
                heuristic_cost(i)=dvmax/sequence_length;
            else
                heuristic_cost(i)=1000;
            
            end
        end
        for dt=100:100:2*T(i);
            k=k+1;
            ma_deg_target=(u^2/h(1,i)^3)*(1-e(1,i)^2)^(1.5)*(t_initial(1,i)+twait+dt)*180/pi;
            ta_target=ma_ta(ma_deg_target,e(1,i));
            dtheta=abs(ta_chaser-ta_target);
            [r2,v2,alpha,universal_anomaly_target]=find_r2_v2(r2_ini,v2_ini,dt,T(i));
            r2=transpose(r2);
            v2=transpose(v2);
            %[dv_prograde_1,dv_retrograde_1]=main_code_book_sub1(r1,r2,dt,dtheta,v1,v2);
            if method==1;
                [v1_prograde,v2_prograde,RAAN_prograde,inclination_prograde,perigee_prograde,true_anomaly_1_prograde,true_anomaly_2_prograde,v1_retrograde,v2_retrograde,RAAN_retrograde,inclination_retrograde,perigee_retrograde,true_anomaly_1_retrograde,true_anomaly_2_retrograde] = lambert_book(r1,r2,dt);
            elseif method==2
                [v1_prograde,v2_prograde,at_short,RAAN_short,inclination_short,perigee_short,v1_retrograde,v2_retrograde,at_long,RAAN_long,inclination_long,perigee_long] = lambert_prussing(r1,r2,dt,dtheta);
            else
                r1=transpose(r1);
                r2=transpose(r2);
                v1=transpose(v1);
                v2=transpose(v2);
                m=1;
                [V1, V2, extremal_distances, exitflag] = lambert(r1, r2, dt, m, u);
                v1_prograde=V1;
                v2_prograde=V2;
                v1_retrograde=V1;
                v2_retrograde=V2;
                
            end
            
            dv1_pro=v1_prograde-v1;
            dv2_pro=v2-v2_prograde;
            dv_prograde_1=norm(dv1_pro)+norm(dv2_pro);
        
            dv1_retro=v1_retrograde-v1;
            dv2_retro=v2-v2_retrograde;
            dv_retrograde_1=norm(dv1_retro)+norm(dv2_retro); 
            
            if isnan(dv_prograde_1)==0;
                dv(i,k)=dv_prograde_1;
                transfer_time(i,k)=dt;
            elseif isnan(dv_retrograde_1)==0; 
                dv(i,k)=dv_retrograde_1;
                transfer_time(i,k)=dt;
            else
                dv(i,k)=nan;
                transfer_time(i,k)=nan;
            end
        end
    else 
        dv(i,1)=1000;
        transfer_time(i,1)=1000;
        heuristic_cost(i)=1000;
    end
end
%dv
%transfer_time
for i=1:tn;
    if i~=sequence_best;
        temp=dv(i,:);
        [minimum]=min(temp(temp>0));
        [x,y]=size(temp);
        for i1=1:y
            if temp(1,i1)==minimum;
                index=i1;
            end
        end
        if isempty(minimum)==0;            
            %i
            %index
            mindv(i)=minimum;
            mindv_dt(i)=transfer_time(i,index);
        else
            mindv(i)=1000;
            mindv_dt(i)=0;
        end
    else
        mindv(i)=1000;
        mindv_dt(i)=1000;
    end
end
%mindv
%mindv_dt
if approach==1;
    dv=[];
    transfer_time=[];
    [dv,index]=min(mindv(mindv>0));
    transfer_time=mindv_dt(index);
    time_elapsed=transfer_time+twait;
    next_node=index;
else
    for i=1:tn;
        totalcost(i)=heuristic_cost(i)+mindv(i);
    end
    [m,index]=min(totalcost);
    dv=[];
    transfer_time=[];
    transfer_time=mindv_dt(index);
    %transfer_time
    time_elapsed=transfer_time+twait;
    dv=mindv(index);
    next_node=index;
end


            
    
    






             
    
    