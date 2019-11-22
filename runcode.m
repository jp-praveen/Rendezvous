function [sequence_1,dv_1,time_elapsed_1,sequence_2,dv_2,time_elapsed_2, sequence_3,dv_3,time_elapsed_3, sequence_4,dv_4,time_elapsed_4]=runcode(dvmax)
sequence_1=[];
sequence_2=[];
sequence_3=[];
sequence_4=[];
dv_1=[];
dv_2=[];
dv_3=[];
dv_4=[];
time_elapsed_1=[];
time_elapsed_2=[];
time_elapsed_3=[];
time_elapsed_4=[];

for j=1:15;
    j
    for i=1:4;
        i
        if i==1;
            [sequence_01,dv_01,time_elapsed_01]=greedy(50,j,dvmax,i)
            sequence_1=padconcatenation(sequence_1,sequence_01,1);
            dv_1=padconcatenation(dv_1,dv_01,1);
            time_elapsed_1=padconcatenation(time_elapsed_1,time_elapsed_01,1);
        elseif i==2;
            [sequence_02,dv_02,time_elapsed_02]=greedy(50,j,dvmax,i)
            sequence_2=padconcatenation(sequence_2,sequence_02,1);
            dv_2=padconcatenation(dv_2,dv_02,1);
            time_elapsed_2=padconcatenation(time_elapsed_2,time_elapsed_02,1);
        elseif i==3;
            [sequence_03,dv_03,time_elapsed_03]=greedy(50,j,dvmax,i)
            sequence_3=padconcatenation(sequence_3,sequence_03,1);
            dv_3=padconcatenation(dv_3,dv_03,1);
            time_elapsed_3=padconcatenation(time_elapsed_3,time_elapsed_03,1);
        else;
            [sequence_04,dv_04,time_elapsed_04]=greedy(50,j,dvmax,i)
            sequence_4=padconcatenation(sequence_4,sequence_04,1);
            dv_4=padconcatenation(dv_4,dv_04,1);
            time_elapsed_4=padconcatenation(time_elapsed_4,time_elapsed_04,1);
        end
    end
end
sequence_1
sequence_2
sequence_3
sequence_4
dv_1
dv_2
dv_3
dv_4
time_elapsed_1
time_elapsed_2
time_elapsed_3
time_elapsed_4
        