%% FIRST BRANCH: multiple CTs
%single cell scenario: graph-coloring
SINRth_c_with_gap = 7 dB; %(>6dB)
SINRth_d_with_gap = 16 dB; %(>15dB)
C = randsample(5,1); %C is random number from 1 to 5 ,4
c=1:C; %CTs, each has its own RB
D=50; %D2D pairs
d=1:D; 
RB(d)=0; %all d are not active initially 
%CTs and DTs are randomly distributed over the cell

d=1;
for c=1:C    
    SINR(c);
    SINR(d);  
    while SINR(c) >= SINRth_c_with_gap
        while 1 % loop is running until some d is assigned to c
            if (SINR(d) >= SINRth_d_with_gap)
                RB(d) = RB(c);
                if (d ==D) %the last pair
                    end program       
                else
                    d++;
                    break while
            else
                RB(d)=0; %turn this DT off   
                if (d ==D) %the last pair
                   end program 
                else
                   d++;
            end
        end while
        SINR(c);
    end while   
end for    
 
%3-cell scenario: run single cell graph-coloring 3 times indepenently,
%run border conditions check C times

%% SECOND BRANCH: mobility
cell_radius=200;
%move all devices randomly at the distance 
            %from 0 to 0.06*cell_radius of their initial place
%whenever any device moved at the distance higher than 0.05*cell_radius - 
            %rerun the whole algorithm, get SINRs, SE
%move all devices randomly at the distance from 0 to 0.04*cell_radius 
            %of their initial place, get SINRs, SE
%compare the performance: how much worse it became?