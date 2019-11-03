% Greedy Approach 
function [sequence_best,dv,time_elapsed] = greedy(totalnode,startnode,dvmax) 
tn=totalnode;
s=startnode;
dv=0;
sequence_best=[s];
time_elapsed=0;
twait=0;
while dv<dvmax;
    [time_elapsed,dv_current,next_node] = realistic_position(tn,s,twait,sequence_best);
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
  