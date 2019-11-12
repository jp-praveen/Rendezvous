function [path] = djkstra(n,edge)               % Changes to be made: 
edge_size=size(edge);                           % if starting and ending nodes are same
m=edge_size(1,1);                               % for any starting and ending node

edge_1=[];
edge_sort=[];
for j=1:n                                        %sorting the edge array in the form of ascendin nodes
    for i=1:m;
        if edge(i,1)==j & edge(i,2)>j ;
            edge_1(:,1)=edge(i,1);
            edge_1(:,2)=edge(i,2);
            edge_1(:,3)=edge(i,3);
            edge_sort=[edge_sort;edge_1];
            
        elseif edge(i,2)==j & edge(i,2)>j;
            edge_1(:,1)=edge(i,2);
            edge_1(:,2)=edge(i,1);
            edge_1(:,3)=edge(i,3);
            edge_sort=[edge_sort;edge_1];
        end
    end
end
edge_sort 
    
 for i=1:n;                                    %Initiating the past costs, zero for node 1 and infinity  for higher nodes
    if i==1;
        past_cost(i,1)=0;
        
    else;
        past_cost(i,1)=500;
    end
 end

    
for i=1:n;
    table(1,i)=past_cost(i,1);                  % Creating the  table with the each column representing each node
    table(2,i)=0;                               % and 1st row corresponds to past cost and 2nd row corresponds to  
end                                             % parent node   
check=[];    
open=[];
close=[];
node_col=edge_sort(:,1);

for j=1:m;                                      % starting the planning from 1st node
    if edge_sort(j,1)==1;
        flag=edge_sort(j,2);
        table(1,flag)=edge_sort(j,3);
        table(2,flag)=edge_sort(j,1);
        check=flag;
        open=[open flag];
    end
end

check1=[];                                      % Making the open array to be in ascending order           
[dum,open_size]=size(open);
for i=1:open_size;
    s=open(1,i);
    check1=[check1 table(1,s)];
end

check1=sort(check1);


for i=1:open_size;
    for j=1:open_size;
        k=open(1,j);                    
        if check1(1,i)==table(1,k);
            check1(1,i)=k;
        end
    end
end
%check1;
open=check1;

current=open(1,1);
f2=1;
while isempty(open)~=1;                                            % Using Dijkstra Algorithm
    for j=1:m;
        if edge_sort(j,1)==current;
            flag=edge_sort(j,2);
            if table(2,flag)==0;
                table(1,flag)=table(1,edge_sort(j,1))+edge_sort(j,3);
                table(2,flag)=edge_sort(j,1);
                check_rep=find(open==flag);
                    if isempty(check_rep)==1;
                        open=[open flag];
                    end
            elseif table(2,flag)~=0;
                if table(1,edge_sort(j,1))+edge_sort(j,3)< table(1,flag);
                    table(1,flag)=table(1,edge_sort(j,1))+edge_sort(j,3);
                    table(2,flag)=edge_sort(j,1);
                    check_rep=find(open==flag);
                    if isempty(check_rep)==1;
                        open=[open flag];
                    end
                else
                    table(1,flag)=table(1,flag);
                end
            end
        end
    end
    close=[close open(1,1)];                           % Including the completed nodes in close array
    [f1,f2]=size(close);
    open(1)=[];                                        % Removing the completed node from open array
    check1=[];
    [dum,open_size]=size(open); 
    for i=1:open_size;                                 % Once again sorting the Open array in ascending order
        s=open(1,i);
        check1=[check1 table(1,s)];
    end
    check1=sort(check1);

    for i=1:open_size;
        for j=1:open_size;
            k=open(1,j);
            if check1(1,i)==table(1,k);
                check1(1,i)=k;
            end
        end
    end
    open=check1;
    if isempty(open)~=1;
        current=open(1,1);
    end
    
end

path=table;

                
             
            
            
      
    


    
    
   
    
            
            
            
    
    