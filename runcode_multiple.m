function [sequence_1,dv_1,time_elapsed_1]=runcode_multiple(dvmax,total_vehicles)
sequence_1=[];
dv_1=[];
time_elapsed_1=[];
previous_sequence=[];

n=total_vehicles;
for i=1:n;
    previous_sequence
    prompt='Enter the start node';
    start_node=input(prompt);
    while ismember(start_node,previous_sequence)==1;
        prompt='Enter the start node';
        start_node=input(prompt);
    end
    [sequence_01,dv_01,time_elapsed_01]=greedy(50,start_node,dvmax,1,previous_sequence)
    previous_sequence=[previous_sequence sequence_01];
    sequence_1=padconcatenation(sequence_1,sequence_01,1);
    dv_1=padconcatenation(dv_1,dv_01,1);
    time_elapsed_1=padconcatenation(time_elapsed_1,time_elapsed_01,1);
end