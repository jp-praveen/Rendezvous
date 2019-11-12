function [opt_path]= a_star(node_xy,obs_xy,edge,opt_dis)
size_node=size(node_xy);
n=size_node(1,1);

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

xn=node_xy(:,2);
yn=node_xy(:,3);
xobs=obs_xy(:,1);
yobs=obs_xy(:,2);
coor_obs=[xobs yobs];
dia_obs=obs_xy(:,3);
rad_obs=dia_obs/2;
scatter(xn,yn,'filled','markerfacecolor',[0 1 0])
a=[1:n]';
b=num2str(a);
c=cellstr(b);
text(xn+0.02,yn+0.02,c);
axis([-1 1 -1 1])
hold on;
viscircles(coor_obs,rad_obs)
for i=1:m;
    n1=edge_sort(i,1);
    n2=edge_sort(i,2);
    x=[xn(n1,1);xn(n2,1)];
    y=[yn(n1,1);yn(n2,1)];
    plot(x,y,'b--');
end  

for i=1:n;                                    %Initiating the past costs, zero for node 1 and infinity  for higher nodes
    if i==1;
        past_cost(i,1)=0;
        
    else;
        past_cost(i,1)=500;
    end
 end
    
for i=1:n;
    table(1,i)=past_cost(i,1);
    table(2,i)=opt_dis(i,1);
    table(3,i)= table(1,i)+ table(2,i);
    table(4,i)=0;
end
       
open=[];
close=[];

for j=1:m;                                      % starting the planning from 1st node
    if edge_sort(j,1)==1;
        flag=edge_sort(j,2);
        table(1,flag)=edge_sort(j,3);
        table(3,flag)= table(1,flag) + table(2,flag);
        table(4,flag)=edge_sort(j,1);
        open=[open flag];
    end
end

check1=[];                                      % Making the open array to be in ascending order           
[dum,open_size]=size(open);
for i=1:open_size;
    s=open(1,i);
    check1=[check1 table(3,s)];
end
check1=sort(check1);
for i=1:open_size;
    for j=1:open_size;
        k=open(1,j);                    
        if check1(1,i)==table(3,k);
            check1(1,i)=k;
        end
    end
end

open=check1;
current=open(1,1);

while isempty(open)~=1;                                            % Using A* Algorithm
    for j=1:m;
        if edge_sort(j,1)==current;
            flag=edge_sort(j,2);
            if table(4,flag)==0;
                table(1,flag)=table(1,edge_sort(j,1))+edge_sort(j,3);
                table(3,flag)=table(1,flag)+table(2,flag);
                table(4,flag)=edge_sort(j,1);
                check_rep=find(open==flag);
                    if isempty(check_rep)==1;
                        open=[open flag];
                    end
            elseif table(2,flag)~=0;
                if table(3,edge_sort(j,1))+edge_sort(j,3)< table(3,flag);
                    table(1,flag)=table(1,edge_sort(j,1))+edge_sort(j,3);
                    table(3,flag)=table(1,flag)+table(2,flag);
                    table(4,flag)=edge_sort(j,1);
                    check_rep=find(open==flag);
                    if isempty(check_rep)==1;
                        open=[open flag];
                    end
                else
                    table(1,flag)=table(1,flag);
                    table(3,flag)=table(1,flag)+table(2,flag);
                end
            end
        end
    end
    close=[close open(1,1)];                           % Including the completed nodes in close array
    open(1)=[];                                        % Removing the completed node from open array
    check1=[];
    [dum,open_size]=size(open); 
    for i=1:open_size;                                 % Once again sorting the Open array in ascending order
        s=open(1,i);
        check1=[check1 table(3,s)];
    end
    check1=sort(check1);

    for i=1:open_size;
        for j=1:open_size;
            k=open(1,j);
            if check1(1,i)==table(3,k);
                check1(1,i)=k;
            end
        end
    end
    open=check1;
    if isempty(open)~=1;
        current=open(1,1);
    end
    
end
opt_path=[];
i=n;
while i~=0;
    opt_path=[opt_path i];
    i=table(4,i);
end

[j,i]=size(opt_path);
for k=1:i-1;
    n1=opt_path(1,k);
    n2=opt_path(1,k+1);
    x=[node_xy(n1,2) node_xy(n2,2)];
    y=[node_xy(n1,3) node_xy(n2,3)];
    plot(x,y,'r')
end


    
    

            
            
            
            
