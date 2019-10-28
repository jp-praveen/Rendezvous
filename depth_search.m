% Depth First Search
function [path] = depth_search(totalnode,startnode,dvmax)        
load matrixmat.mat;
matrixmat=matrixmat{:,:};
matrixmat(1,:) = [];
matrixmat(:,1) = [];

tn=totalnode;
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
    if temp(1)==matrixmat(s1,i);
        open1=[open1 i];
        temp(1)=[];
        if isempty(temp)==1;
            break;
        end
    end
end
open1
%cost1
[m,n]=size(open1);

while isempty(open1)==0;
    s2=open1(1);
    open2=[];
    close1=[];
    close2=[];
    close3=[];
    close4=[];
    open3=[];
    open4=[];
    open5=[];
    sequence=[];
    total_cost=[];
    k=3;
    table_info=[];
    for j=1:tn;
        if isnan(matrixmat(s2,j))==0 & cost1(1)+matrixmat(s2,j)<dvmax;
            open2=[open2 j];
            close1=[close1 open1(1)];
            table_info(1,j)=cost1(1)+matrixmat(s2,j);
            table_info(2,j)=s2;
            table_info(3,j)=k;
            temp_seq=[s1 open1(1) j 0 0 0 0 0 0 0 ];
            sequence= vertcat(sequence,temp_seq);
            total_cost=vertcat(total_cost,table_info(1,j));
            for q=1:tn;
                s3=max(open2(1));
                close2=[close2 open2(1)];
                %open2=[];
                if isnan(matrixmat(s3,q))==0 & table_info(1,j)+matrixmat(s3,q)<dvmax;
                    open3=[open3 q];
                    table_info(1,q)=table_info(1,j)+matrixmat(s3,q);
                    table_info(2,q)=s3;
                    table_info(3,q)=k+1;
                    temp_seq=[s1 open1(1) j q 0 0 0 0 0 0];
                    sequence= vertcat(sequence,temp_seq);
                    total_cost=vertcat(total_cost,table_info(1,q));
                    for w=1:tn;
                        s4=max(open3(1));
                        close3=[close3 open3(1)];
                        %open2=[];
                        if isnan(matrixmat(s4,w))==0 & table_info(1,q)+matrixmat(s4,w)<dvmax;
                            open4=[open4 w];
                            table_info(1,w)=table_info(1,q)+matrixmat(s4,w);
                            table_info(2,w)=s4;
                            table_info(3,w)=k+2;
                            temp_seq=[s1 open1(1) j q w 0 0 0 0 0];
                            sequence= vertcat(sequence,temp_seq);
                            total_cost=vertcat(total_cost,table_info(1,w));
                            for e=1:tn;
                                s5=max(open4(1));
                                close4=[close4 open4(1)];
                                %open3=[];
                                if isnan(matrixmat(s5,e))==0 & table_info(1,w)+matrixmat(s5,e)<dvmax;
                                    open5=[open5 w];
                                    table_info(1,e)=table_info(1,w)+matrixmat(s5,e);
                                    table_info(2,e)=s5;
                                    table_info(3,e)=k+3;     
                                    temp_seq=[s1 open1(1) j q w e 0 0 0 0 ];
                                    sequence= vertcat(sequence,temp_seq);
                                    total_cost=vertcat(total_cost,table_info(1,e));
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
     sequence
     total_cost
end

                    


        


        