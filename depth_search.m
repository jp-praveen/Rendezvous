% Depth First Search
function [sequence_best] = depth_search(totalnode,startnode,dvmax)        

%function [possible_start,possible_end]=sf1_depth_search(totalnode)
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
%end
s1=startnode;
open1=[];
cost1=[];
for i=1:tn;
    if isnan(matrixmat(s1,i))==0;
        open1=[open1 i];
        cost1=[cost1 matrixmat(s1,i)];
    end
end
%open1
%cost1
[m,n]=size(open1);
I=find(cost1<dvmax);
temp=[];
while I~=0;
    temp=[temp cost1(I(1))];
    I(1)=[];
end
cost1=temp;
open1=[];
for i=1:tn;
    if isempty(temp)==1;
        break;
    end
    if temp(1)==matrixmat(s1,i);
        open1=[open1 i];
        temp(1)=[];
        if isempty(temp)==1;
            break;
        end
    end
end
%open1
%cost1
[m,n]=size(open1);
sequence_all=[];
total_cost_all=[];
while isempty(open1)==0;
    s2=open1(1);
    open2=[];
    close1=[];
    close2=[];
    close3=[];
    close4=[];
    close5=[];
    close6=[];
    close7=[];
    open3=[];
    open4=[];
    open5=[];
    open6=[];
    open7=[];
    open8=[];
    sequence=[];
    total_cost=[];
    k=3;
    table_info=[];
    for j=1:tn;
        if isnan(matrixmat(s2,j))==0 & cost1(1)+matrixmat(s2,j)<dvmax &j~=s2 &j~=s1;
            open2=[open2 j];
            close1=[close1 s2];
            table_info(1,j)=cost1(1)+matrixmat(s2,j);
            table_info(2,j)=s2;
            table_info(3,j)=k;
            temp_seq=[s1 open1(1) j 0 0 0 0 0 0 0 ];
            sequence= vertcat(sequence,temp_seq);
            total_cost=vertcat(total_cost,table_info(1,j));
            for q=1:tn;
                s3=max(open2);
                close2=[close2 s3];
                %open2=[];
                if isnan(matrixmat(s3,q))==0 & table_info(1,j)+matrixmat(s3,q)<dvmax & q~=s1 & q~=open1(1) &q~=j  ;
                    open3=[open3 q];
                    table_info(1,q)=table_info(1,j)+matrixmat(s3,q);
                    table_info(2,q)=s3;
                    table_info(3,q)=k+1;
                    temp_seq=[s1 open1(1) j q 0 0 0 0 0 0];
                    sequence= vertcat(sequence,temp_seq);
                    total_cost=vertcat(total_cost,table_info(1,q));
                    for w=1:tn;
                        s4=max(open3);
                        close3=[close3 s4];
                        %open2=[];
                        if isnan(matrixmat(s4,w))==0 & table_info(1,q)+matrixmat(s4,w)<dvmax & w~=s1 & w~=open1(1) & w~=j &w~=q;
                            open4=[open4 w];
                            table_info(1,w)=table_info(1,q)+matrixmat(s4,w);
                            table_info(2,w)=s4;
                            table_info(3,w)=k+2;
                            temp_seq=[s1 open1(1) j q w 0 0 0 0 0];
                            sequence= vertcat(sequence,temp_seq);
                            total_cost=vertcat(total_cost,table_info(1,w));
                            for e=1:tn;
                                s5=max(open4);
                                close4=[close4 s5];
                                %open3=[];
                                if isnan(matrixmat(s5,e))==0 & table_info(1,w)+matrixmat(s5,e)<dvmax & e~=s1 & e~=open1(1) & e~=j & e~=q &e~=w ;
                                    open5=[open5 w];
                                    table_info(1,e)=table_info(1,w)+matrixmat(s5,e);
                                    table_info(2,e)=s5;
                                    table_info(3,e)=k+3;     
                                    temp_seq=[s1 open1(1) j q w e 0 0 0 0 ];
                                    sequence= vertcat(sequence,temp_seq);
                                    total_cost=vertcat(total_cost,table_info(1,e));
                                    for r=1:tn;
                                        s6=max(open5);
                                        close5=[close5 s6];
                                        if isnan(matrixmat(s6,r))==0 & table_info(1,e)+matrixmat(s6,r)<dvmax & r~=s1 & r~=open1(1) & r~=j & r~=q &r~=w &r~=e;
                                            open6=[open6 e];
                                            table_info(1,r)=table_info(1,e)+matrixmat(s6,r);
                                            table_info(2,r)=s6;
                                            table_info(3,r)=k+4;     
                                            temp_seq=[s1 open1(1) j q w e r 0 0 0 ];
                                            sequence= vertcat(sequence,temp_seq);
                                            total_cost=vertcat(total_cost,table_info(1,r));
                                            for t=1:tn;
                                                s7=max(open6);
                                                close6=[close6 s7];
                                                if isnan(matrixmat(s7,t))==0 & table_info(1,r)+matrixmat(s7,t)<dvmax & t~=s1 & t~=open1(1) & t~=j & t~=q &t~=w &t~=e &t~=r;
                                                    open7=[open7 r];
                                                    table_info(1,t)=table_info(1,r)+matrixmat(s7,t);
                                                    table_info(2,t)=s7;
                                                    table_info(3,t)=k+5;     
                                                    temp_seq=[s1 open1(1) j q w e r t 0 0 ];
                                                    sequence= vertcat(sequence,temp_seq);
                                                    total_cost=vertcat(total_cost,table_info(1,t));
                                                    for y=1:tn;
                                                        s8=max(open7);
                                                        close7=[close7 s8];
                                                        if isnan(matrixmat(s8,y))==0 & table_info(1,t)+matrixmat(s8,y)<dvmax & y~=s1 & y~=open1(1) & y~=j & y~=q &y~=w &y~=e &y~=r &y~=t;
                                                            open8=[open8 t];
                                                            table_info(1,y)=table_info(1,t)+matrixmat(s8,y);
                                                            table_info(2,y)=s8;
                                                            table_info(3,y)=k+6;     
                                                            temp_seq=[s1 open1(1) j q w e r t y 0 ];
                                                            sequence= vertcat(sequence,temp_seq);
                                                            total_cost=vertcat(total_cost,table_info(1,y));
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                                          
                                end
                            end
                        end
                    end
                end
            end
        end
        %open4=[];
    end
    cost1(1)=[];
    open1(1)=[];
     table_info;
     sequence;
     sequence_all=[sequence_all;sequence];
     total_cost_all=[total_cost_all;total_cost]; 
     
     total_cost;
     %sequence_cost=[sequence total_cost]
     %sequence_best=[sequence_best cost_best]
     
         
end
sequence=sequence_all;
total_cost=total_cost_all;
sequence_best=[];
cost_best=[];
[m,n]=size(sequence);
for k=1:m;
    temp=[];
    for l=1:10
        if sequence(k,l)~=0;
            temp=[temp sequence(k,l)];
        else;
            break;
        end
    end
    if length(temp)>length(sequence_best);
        sequence_best=temp;
        cost_best=total_cost(k);
    elseif length(temp)==length(sequence_best);
        if cost_best>total_cost(k);
            sequence_best=temp;
            cost_best=total_cost(k);
        end
    end
end
sequence_best=[sequence_best cost_best];
end



        


        