clc;clear;close all;format compact;

CUE_SINR_min_dB = 20; % CT SINR th, level of forb reg harshness, minor influence
DUE_SINR_min_dB	= 40; %DT SINR th     6/15 ARP, 30/46 opt
CUE_SINR_min = db2pow(CUE_SINR_min_dB); %times
DUE_SINR_min = db2pow(DUE_SINR_min_dB); %times
CUE_Exp = 3.67; % Cellular Pathloss Exponent - Alpha_zero
DUE_Exp = 4.33; % Devices Pathloss Exponent - Alpha_d

%coordinates of BSs, 3 cells
eNB1_x = 200; eNB1_y = 200; 
eNB2_x = 400; eNB2_y = 200+200*sqrt(3);
eNB3_x = 600; eNB3_y = 200;
customColormap = [1 0 0  ;1 0 1 ;0 0 1; 1 1 0; 0 1 1 ; 0 1 0]; %r m c y b g
customColormap_char = ['r'; 'm'; 'b' ;'y'; 'c'; 'g' ]; 

mobility=5; % percent of the cell radius
velocity_kmh=10; %km/h
velocity=velocity_kmh*5/18;
timestamp=10; %sec
number_of_timestamps=2;
CTs = 3; % same # of CTs for each cell
no_runs=200;

for ct=1:CTs
    SINR_C_mAOS{1,ct}=[];
    SINR_D_mAOS{1,ct}=[];
    SINR_C_mAOS{2,ct}=[];
    SINR_D_mAOS{2,ct}=[];
    SINR_C_mAOS_dB{1,ct}=[];
    SINR_D_mAOS_dB{1,ct}=[];
    SINR_C_mAOS_dB{2,ct}=[];
    SINR_D_mAOS_dB{2,ct}=[];
    CT_BS_gain1{1,ct}=[];CT_BS_gain2{1,ct}=[];CT_BS_gain3{1,ct}=[];
    CT_BS_gain1{2,ct}=[];CT_BS_gain2{2,ct}=[];CT_BS_gain3{2,ct}=[];
    SINR_C_rand_all{1,ct}=[];SINR_C_rand_all{2,ct}=[];
    SINR_D_rand_all{1,ct}=[];SINR_D_rand_all{2,ct}=[];
    SINR_C_rand_all_dB{1,ct}=[];SINR_C_rand_all_dB{2,ct}=[];
    SINR_D_rand_all_dB{1,ct}=[];SINR_D_rand_all_dB{2,ct}=[];
end

Cell_Radius = 200;
D2D_Sep_Max = 0.1*Cell_Radius;
Max_Users = 100; % pairs

for n=1:no_runs
    to_disp=['n# ',num2str(n)];
    disp(to_disp);
    Max_Users1=Max_Users; Max_Users2=Max_Users; Max_Users3=Max_Users;
    
    %distibution of users in 3 cells
    D2D_user_list1 = LTE_UE_uniform_distribution_upd(eNB1_x,eNB1_y,Cell_Radius,D2D_Sep_Max, Max_Users); %Tx_x,Tx_y,Rx_x,Rx_y
    D2D_user_list2 = LTE_UE_uniform_distribution_upd(eNB2_x,eNB2_y,Cell_Radius,D2D_Sep_Max, Max_Users);
    D2D_user_list3 = LTE_UE_uniform_distribution_upd(eNB3_x,eNB3_y,Cell_Radius,D2D_Sep_Max, Max_Users);
    [CUEs{1,1},~, ~,~, ~,~, ~,...
        CT_BS_gain1{1,1},CT_BS_gain2{1,1},CT_BS_gain3{1,1}]=generate_3_CTs(eNB1_x,eNB1_y,eNB2_x,eNB2_y,eNB3_x,eNB3_y,Cell_Radius);
    
    D2D_user_list1_initial=D2D_user_list1;
    D2D_user_list2_initial=D2D_user_list2;
    D2D_user_list3_initial=D2D_user_list3;
    
    [D2D_user_list1_moved]=move_list(D2D_user_list1_initial,mobility,Cell_Radius); %move DTs
    [D2D_user_list2_moved]=move_list(D2D_user_list2_initial,mobility,Cell_Radius);
    [D2D_user_list3_moved]=move_list(D2D_user_list3_initial,mobility,Cell_Radius);

          
    %N CTs
    for ct=1:CTs
        switch ct
            case 1
                color=customColormap_char(1); %r
            case 2
                color=customColormap_char(2); %m
            case 3
                color=customColormap_char(3); %c
        end
        
        %RANDOM
        rand_sqc1 = randperm(size(D2D_user_list1,1));  %https://www.mathworks.com/matlabcentral/answers/17659-how-to-select-random-rows-from-a-matrix
        random_list1{1,ct} = D2D_user_list1(rand_sqc1(1:5),:);
        rand_sqc2 = randperm(size(D2D_user_list2,1));
        random_list2{1,ct} = D2D_user_list2(rand_sqc2(1:5),:);
        rand_sqc3 = randperm(size(D2D_user_list3,1));
        random_list3{1,ct} = D2D_user_list3(rand_sqc3(1:5),:);
        random_list_all{1,ct}=[random_list1{1,ct};random_list2{1,ct};random_list3{1,ct}];
        pairs_qtity_random{1,ct}(n,1)=size(random_list1{1,ct},1)+size(random_list2{1,ct},1)+size(random_list3{1,ct},1); 
    
        [SINR_C1_rand{1,ct},SINR_D1_rand{1,ct}] = SINR_calc (1,random_list1{1,ct},eNB1_x,eNB1_y,CUEs{1,ct}(1,1),CUEs{1,ct}(1,2),CUEs{1,ct}(2,1),CUEs{1,ct}(2,2),CUEs{1,ct}(3,1),CUEs{1,ct}(3,2),random_list_all{1,ct},CT_BS_gain1{1,ct},CT_BS_gain2{1,ct},CT_BS_gain3{1,ct});
        [SINR_C2_rand{1,ct},SINR_D2_rand{1,ct}] = SINR_calc (2,random_list2{1,ct},eNB2_x,eNB2_y,CUEs{1,ct}(1,1),CUEs{1,ct}(1,2),CUEs{1,ct}(2,1),CUEs{1,ct}(2,2),CUEs{1,ct}(3,1),CUEs{1,ct}(3,2),random_list_all{1,ct},CT_BS_gain1{1,ct},CT_BS_gain2{1,ct},CT_BS_gain3{1,ct});
        [SINR_C3_rand{1,ct},SINR_D3_rand{1,ct}] = SINR_calc (3,random_list3{1,ct},eNB3_x,eNB3_y,CUEs{1,ct}(1,1),CUEs{1,ct}(1,2),CUEs{1,ct}(2,1),CUEs{1,ct}(2,2),CUEs{1,ct}(3,1),CUEs{1,ct}(3,2),random_list_all{1,ct},CT_BS_gain1{1,ct},CT_BS_gain2{1,ct},CT_BS_gain3{1,ct});
        mult=1;
        for k=1:size(SINR_D1_rand{1,ct},1)
            mult=mult * ( 1 + SINR_D1_rand{1,ct}(k) );
        end
        for k=1:size(SINR_D2_rand{1,ct},1)
            mult=mult * ( 1 + SINR_D2_rand{1,ct}(k) );
        end
        for k=1:size(SINR_D3_rand{1,ct},1)
            mult=mult * ( 1 + SINR_D3_rand{1,ct}(k) );
        end
        SE_rand{1,ct}(n,:)=log2( (1+SINR_C1_rand{1,ct})*(1+SINR_C2_rand{1,ct})*(1+SINR_C3_rand{1,ct}) * mult );
        SINR_C_rand_all{1,ct}=[SINR_C_rand_all{1,ct};SINR_C1_rand{1,ct};SINR_C2_rand{1,ct};SINR_C3_rand{1,ct}];
        SINR_D_rand_all{1,ct}=[SINR_D_rand_all{1,ct};SINR_D1_rand{1,ct};SINR_D2_rand{1,ct};SINR_D3_rand{1,ct}];
        SINR_C_rand_all_dB{1,ct}=pow2db(SINR_C_rand_all{1,ct}); %random
        SINR_D_rand_all_dB{1,ct}=pow2db(SINR_D_rand_all{1,ct});

        %mAOS
        [AOS_user_list1_upd,AOS_user_list2_upd,AOS_user_list3_upd,...
        D2Ds_left1,D2Ds_left2,D2Ds_left3,...
        AOS_user_list1,AOS_user_list2,AOS_user_list3,...
        SE_new_AOS,SINR_C_new_all_AOS_dB,SINR_D_i_new_all_AOS_dB,numbers_of_pairs_new_AOS,...
        rank_AOS1{1,ct},rank_AOS2{1,ct},rank_AOS3{1,ct}] = ...
        multi_cell_AOS ...
        (D2D_user_list1,D2D_user_list2,D2D_user_list3, Cell_Radius,D2D_Sep_Max,...
        Max_Users1,Max_Users2,Max_Users3, ...
        CUE_Exp,DUE_Exp,CUE_SINR_min,DUE_SINR_min,eNB1_x,eNB2_x,eNB3_x,eNB1_y,eNB2_y,eNB3_y,...
        CUEs{1,ct}(1,1),CUEs{1,ct}(1,2),CUEs{1,ct}(2,1),CUEs{1,ct}(2,2),CUEs{1,ct}(3,1),CUEs{1,ct}(3,2),...
        CT_BS_gain1{1,ct},CT_BS_gain2{1,ct},CT_BS_gain3{1,ct}); 
        
        AOS_list1{1,ct}=AOS_user_list1;
        AOS_list2{1,ct}=AOS_user_list2;
        AOS_list3{1,ct}=AOS_user_list3;
        mAOS_list1{1,ct}=AOS_user_list1_upd; %save selected mAOS pairs
        mAOS_list2{1,ct}=AOS_user_list2_upd;
        mAOS_list3{1,ct}=AOS_user_list3_upd;
        mAOS_list_all{1,ct}=[mAOS_list1{1,ct};mAOS_list2{1,ct};mAOS_list3{1,ct}]; %each col here consists of all DTs in 3 cells sharing same RB (1 to 3)

        if (n==1)
            distr_visual(color,CUEs{1,ct}(1,1),CUEs{1,ct}(1,2),CUEs{1,ct}(2,1),CUEs{1,ct}(2,2),CUEs{1,ct}(3,1),CUEs{1,ct}(3,2),... %plot the distribution
            D2D_user_list1,D2D_user_list2,D2D_user_list3,...
        random_list1{1,ct},random_list2{1,ct},random_list3{1,ct},...
            mAOS_list1{1,ct},mAOS_list2{1,ct},mAOS_list3{1,ct});
           %AOS_list1{1,ct},AOS_list2{1,ct},AOS_list3{1,ct},...
        end
        
        SE_mAOS{1,ct}(n,1)=SE_new_AOS;
        pairs_qtity_mAOS{1,ct}(n,1)=numbers_of_pairs_new_AOS;
        SINR_C_mAOS_dB{1,ct}=[SINR_C_mAOS_dB{1,ct};SINR_C_new_all_AOS_dB];
        SINR_D_mAOS_dB{1,ct}=[SINR_D_mAOS_dB{1,ct};SINR_D_i_new_all_AOS_dB];

             
        rank_mAOS1{1,ct}= ismember(D2D_user_list1_initial,mAOS_list1{1,ct}); %find indeces of the selected pairs
        rank_mAOS2{1,ct}= ismember(D2D_user_list2_initial,mAOS_list2{1,ct}); %after moving this selection matrix remains the same
        rank_mAOS3{1,ct}= ismember(D2D_user_list3_initial,mAOS_list3{1,ct});
        
       
     
    %mobility: 1 hop
        for cell=1:3 %move CT in each of 3 cells 
            theta_= 2*pi*rand(1,1);
            moving_dist=mobility*0.01*Cell_Radius*sqrt(rand(1,1));
            CUEs{2,ct}(cell,1) = moving_dist*cos(theta_) + CUEs{1,ct}(cell,1) ;%x
            CUEs{2,ct}(cell,2) = moving_dist*sin(theta_) + CUEs{1,ct}(cell,2) ;%y
        end
        CT_BS_gain1{2,ct} = LTE_channel_model_urban_micro_NLOS(CUEs{2,ct}(1,1),CUEs{2,ct}(1,2),eNB1_x,eNB1_y); 
        CT_BS_gain2{2,ct} = LTE_channel_model_urban_micro_NLOS(CUEs{2,ct}(2,1),CUEs{2,ct}(2,2),eNB2_x,eNB2_y);
        CT_BS_gain3{2,ct} = LTE_channel_model_urban_micro_NLOS(CUEs{2,ct}(3,1),CUEs{2,ct}(3,2),eNB3_x,eNB3_y);


        mAOS_list1{2,ct} = mask (D2D_user_list1_moved,rank_mAOS1{1,ct}(:,1)); % select same AOS pairs
        mAOS_list2{2,ct} = mask (D2D_user_list2_moved,rank_mAOS2{1,ct}(:,1)); 
        mAOS_list3{2,ct} = mask (D2D_user_list3_moved,rank_mAOS3{1,ct}(:,1)); 
        mAOS_list_all{2,ct}=[mAOS_list1{2,ct};mAOS_list2{2,ct};mAOS_list3{2,ct}]; 
        
%         test=~(ismember(D2D_user_list1_moved,mAOS_list1{2,ct}));
%         D2D_user_list1_moved_decreased=mask(D2D_user_list1_moved,test(:,1));
%         test=~(ismember(D2D_user_list2_moved,mAOS_list2{2,ct}));
%         D2D_user_list2_moved_decreased=mask(D2D_user_list2_moved,test(:,1));
%         test=~(ismember(D2D_user_list3_moved,mAOS_list3{2,ct}));
%         D2D_user_list3_moved_decreased=mask(D2D_user_list3_moved,test(:,1));
           
        AOS_list1{2,ct} = mask (D2D_user_list1_moved,rank_AOS1{1,ct}); % select same AOS pairs
        AOS_list2{2,ct} = mask (D2D_user_list2_moved,rank_AOS2{1,ct});
        AOS_list3{2,ct} = mask (D2D_user_list3_moved,rank_AOS3{1,ct});

%          if ((ct==1||ct==2)&&n==1) %plot moved distribution visualization          
%             distr_visual(color,CUEs{2,ct}(1,1),CUEs{2,ct}(1,2),CUEs{2,ct}(2,1),CUEs{2,ct}(2,2),CUEs{2,ct}(3,1),CUEs{2,ct}(3,2),... %plot the distribution after moving
%             D2D_user_list1_moved,D2D_user_list2_moved,D2D_user_list3_moved,...
%             [],[],[],...
%            mAOS_list1{2,ct},mAOS_list2{2,ct},mAOS_list3{2,ct});
          %AOS_list1{2,ct},AOS_list2{2,ct},AOS_list3{2,ct},...
%          end

        
        [SINR_C1_mAOS{2,ct},SINR_D1_mAOS{2,ct}] = SINR_calc (1,mAOS_list1{2,ct},eNB1_x,eNB1_y,CUEs{2,ct}(1,1),CUEs{2,ct}(1,2),CUEs{2,ct}(2,1),CUEs{2,ct}(2,2),CUEs{2,ct}(3,1),CUEs{2,ct}(3,2),mAOS_list_all{2,ct},CT_BS_gain1{2,ct},CT_BS_gain2{2,ct},CT_BS_gain3{2,ct});
        [SINR_C2_mAOS{2,ct},SINR_D2_mAOS{2,ct}] = SINR_calc (2,mAOS_list2{2,ct},eNB2_x,eNB2_y,CUEs{2,ct}(1,1),CUEs{2,ct}(1,2),CUEs{2,ct}(2,1),CUEs{2,ct}(2,2),CUEs{2,ct}(3,1),CUEs{2,ct}(3,2),mAOS_list_all{2,ct},CT_BS_gain1{2,ct},CT_BS_gain2{2,ct},CT_BS_gain3{2,ct});
        [SINR_C3_mAOS{2,ct},SINR_D3_mAOS{2,ct}] = SINR_calc (3,mAOS_list3{2,ct},eNB3_x,eNB3_y,CUEs{2,ct}(1,1),CUEs{2,ct}(1,2),CUEs{2,ct}(2,1),CUEs{2,ct}(2,2),CUEs{2,ct}(3,1),CUEs{2,ct}(3,2),mAOS_list_all{2,ct},CT_BS_gain1{2,ct},CT_BS_gain2{2,ct},CT_BS_gain3{2,ct});
            mult=1;
            for k=1:size(SINR_D1_mAOS{2,ct},1)
                mult=mult * ( 1 + SINR_D1_mAOS{2,ct}(k) );
            end
            for k=1:size(SINR_D2_mAOS{2,ct},1)
                mult=mult * ( 1 + SINR_D2_mAOS{2,ct}(k) );
            end
            for k=1:size(SINR_D3_mAOS{2,ct},1)
                mult=mult * ( 1 + SINR_D3_mAOS{2,ct}(k) );
            end
            SE_mAOS{2,ct}(n,:)=log2( (1+SINR_C1_mAOS{2,ct})*(1+SINR_C2_mAOS{2,ct})*(1+SINR_C3_mAOS{2,ct}) * mult );
            %pairs_qtity_mAOS is the same throughout all hops
            SINR_C_mAOS{2,ct}=[SINR_C_mAOS{2,ct};SINR_C1_mAOS{2,ct};SINR_C2_mAOS{2,ct};SINR_C3_mAOS{2,ct}];
            SINR_D_mAOS{2,ct}=[SINR_D_mAOS{2,ct};SINR_D1_mAOS{2,ct};SINR_D2_mAOS{2,ct};SINR_D3_mAOS{2,ct}];
            SINR_C_mAOS_dB{2,ct}=pow2db(SINR_C_mAOS{2,ct});
            SINR_D_mAOS_dB{2,ct}=pow2db(SINR_D_mAOS{2,ct});
            
            
        D2D_user_list1=D2Ds_left1; %!!! should be at the end
        D2D_user_list2=D2Ds_left2; %output will go on the input of the mAOS
        D2D_user_list3=D2Ds_left3;
        
        if ~(ct==CTs) %generate new CTs for new RB if it is not the last CT set for reuse 
            [CUEs{1,ct+1},~, ~,~, ~,~, ~,...  
            CT_BS_gain1{1,ct+1},CT_BS_gain2{1,ct+1},CT_BS_gain3{1,ct+1}]=generate_3_CTs...
            (eNB1_x,eNB1_y,eNB2_x,eNB2_y,eNB3_x,eNB3_y,Cell_Radius);

            Max_Users1=size(D2Ds_left1,1);
            Max_Users2=size(D2Ds_left2,1);
            Max_Users3=size(D2Ds_left3,1);
        end

    end %end CTs 


end %end no_runs     
    


%% PLOTS
%CDF DT
figure
xlim([0 50]);
hold on
for ct=1:CTs
    dotted_SINR_D(ct)=cdfplot(SINR_D_rand_all_dB{1,ct}); %random
    set( dotted_SINR_D(ct), 'LineStyle', ':','Color',customColormap(ct,:));
    
    solid_SINR_D(ct)=cdfplot(SINR_D_mAOS_dB{1,ct}); %mAOS
    set( solid_SINR_D(ct),'Color',customColormap(ct,:));

%     dashed_SINR_D(ct)=cdfplot(SINR_D_mAOS_dB{2,ct}); %mAOS+hop
%     set( dashed_SINR_D(ct), 'LineStyle', '--','Color',customColormap(ct,:));
end
xline(15,'g'); %outage probability
grid on
legend('random RB1','mAOS RB1','random RB2','mAOS RB2','random RB3','mAOS RB3','outage probability');
% legend('random RB1','mAOS RB1','mAOS RB1 + hop','random RB2','mAOS RB2','mAOS RB2 + hop','random RB3','mAOS RB3','mAOS RB3 + hop','outage probability');
xlabel('DT SINR (dB)','FontName','Arial','FontSize',14);
ylabel('CDF','FontName','Arial','FontSize',14);

%CDF CT
figure
xlim([0 50]);
hold on
for ct=1:CTs
    dotted_SINR_C(ct)=cdfplot(SINR_C_rand_all_dB{1,ct}); %random
    set( dotted_SINR_C(ct), 'LineStyle', ':','Color',customColormap(ct,:));

    solid_SINR_C(ct)=cdfplot(SINR_C_mAOS_dB{1,ct});
    set( solid_SINR_C(ct),'Color',customColormap(ct,:));

%     dashed_SINR_C(ct)=cdfplot(SINR_C_mAOS_dB{2,ct});
%     set( dashed_SINR_C(ct), 'LineStyle', '--','Color',customColormap(ct,:));
end
xline(6,'g'); %outage probability
grid on
legend('random RB1','mAOS RB1','random RB2','mAOS RB2','random RB3','mAOS RB3','outage probability');
xlabel('CT SINR (dB)','FontName','Arial','FontSize',14);
ylabel('CDF','FontName','Arial','FontSize',14);

%SE
figure
for ct=1:CTs %append 0,0 value for better representation
   mean_values_random{1,ct}=find_mean_for_each_pairs_qtity(pairs_qtity_random{1,ct},SE_rand{1,ct});
   mean_values_random{1,ct}=[0 0; mean_values_random{1,ct}];
   plot(mean_values_random{1,ct}(:,1),mean_values_random{1,ct}(:,2),':o','linewidth',1,'Color',customColormap(ct,:));
   hold on
   mean_values{1,ct}=find_mean_for_each_pairs_qtity(pairs_qtity_mAOS{1,ct},SE_mAOS{1,ct});
   mean_values{1,ct}=[0 0; mean_values{1,ct}];
   plot(mean_values{1,ct}(:,1),mean_values{1,ct}(:,2),'-o','linewidth',1,'Color',customColormap(ct,:));
   
%    mean_values{2,ct}=find_mean_for_each_pairs_qtity(pairs_qtity_mAOS{1,ct},SE_mAOS{2,ct});
%    mean_values{2,ct}=[0 0; mean_values{2,ct}];
%    plot(mean_values{2,ct}(:,1),mean_values{2,ct}(:,2),'--o','linewidth',1,'Color',customColormap(ct,:));

end
   yline(147,'g'); %av SE 
   xline(16,'g'); %av # of pairs

grid on
legend('random RB1','mAOS RB1','random RB2','mAOS RB2','random RB3','mAOS RB3','av. SE & # of pairs');
% legend('random RB1','mAOS RB1','mAOS RB1 + hop','random RB2','mAOS RB2','mAOS RB2 + hop','random RB3','mAOS RB3','mAOS RB3 + hop','av. SE & # of pairs');
xlabel('Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);





%% FUNCTIONS
function [list_moved]=move_list(list,mobility,Cell_Radius) 
    for i=1:size(list,1)
        theta_= 2*pi*rand(1,1);
        moving_dist=mobility*0.01*Cell_Radius*sqrt(rand(1,1));
        list_moved(i,1) = moving_dist*cos(theta_) + list(i,1) ;
        list_moved(i,2) = moving_dist*sin(theta_) + list(i,2) ;
        theta_= 2*pi*rand(1,1);
        moving_dist=mobility*0.01*Cell_Radius*sqrt(rand(1,1));
        list_moved(i,3) = moving_dist*cos(theta_) + list(i,3) ;
        list_moved(i,4) = moving_dist*sin(theta_) + list(i,4) ;
    end
end
function [values]=find_mean_for_each_pairs_qtity(numbers_of_pairs,SE)
    A=[numbers_of_pairs,SE]; 
    A=sortrows(A); %from smallest # of pairs to highest
    [Alink,aa] = findgroups(A(:,1)); %aa - different values of pairs, Alink - links from aa to A
    values = [ aa, splitapply(@mean,A(:,2),Alink)]; %find mean value for each # of pairs
end
function distr_visual(color, CUE1_x, CUE1_y,CUE2_x, CUE2_y,CUE3_x, CUE3_y,...
    D2D_user_list1,D2D_user_list2,D2D_user_list3,...
    AOS_user_list1,AOS_user_list2,AOS_user_list3,...
    AOS_user_list1_upd,AOS_user_list2_upd,AOS_user_list3_upd)

    figure
    viscircles([200 200],200,'Color','k','linewidth',0.5);
    viscircles([600 200],200,'Color','k','linewidth',0.5);
    viscircles([400 200+200*sqrt(3)],200,'Color','k','linewidth',0.5);
    axis ([-10 810 -10 800]);
    axis square
    grid on;
    hold on;
    scatter([CUE1_x CUE2_x CUE3_x],[CUE1_y CUE2_y CUE3_y],'green','filled','MarkerEdgeAlpha',.2);
    scatter([CUE1_x CUE2_x CUE3_x],[CUE1_y CUE2_y CUE3_y],color);
    
    scatter(D2D_user_list1(:,1),D2D_user_list1(:,2),'k','MarkerEdgeAlpha',.2);
    scatter(D2D_user_list1(:,3),D2D_user_list1(:,4),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
    scatter(D2D_user_list2(:,1),D2D_user_list2(:,2),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
    scatter(D2D_user_list2(:,3),D2D_user_list2(:,4),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
    scatter(D2D_user_list3(:,1),D2D_user_list3(:,2),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
    scatter(D2D_user_list3(:,3),D2D_user_list3(:,4),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
    
    if ~(isempty(AOS_user_list1))
    scatter(AOS_user_list1(:,1),AOS_user_list1(:,2),color);
    scatter(AOS_user_list1(:,3),AOS_user_list1(:,4),color,'HandleVisibility','off');
    end
    if ~(isempty(AOS_user_list2))
    scatter(AOS_user_list2(:,1),AOS_user_list2(:,2),color,'HandleVisibility','off');
    scatter(AOS_user_list2(:,3),AOS_user_list2(:,4),color,'HandleVisibility','off');
    end
    if ~(isempty(AOS_user_list3))
    scatter(AOS_user_list3(:,1),AOS_user_list3(:,2),color,'HandleVisibility','off');
    scatter(AOS_user_list3(:,3),AOS_user_list3(:,4),color,'HandleVisibility','off');
    end
    
    if ~(isempty(AOS_user_list1_upd))
    scatter(AOS_user_list1_upd(:,1),AOS_user_list1_upd(:,2),color,'filled');
    scatter(AOS_user_list1_upd(:,3),AOS_user_list1_upd(:,4),color,'filled','HandleVisibility','off');
    end
    if ~(isempty(AOS_user_list2_upd))
    scatter(AOS_user_list2_upd(:,1),AOS_user_list2_upd(:,2),color,'filled','HandleVisibility','off'); %HandleVisibility is off for legend to have only 3 opaque properties
    scatter(AOS_user_list2_upd(:,3),AOS_user_list2_upd(:,4),color,'filled','HandleVisibility','off');
    end
    if ~(isempty(AOS_user_list3_upd))
    scatter(AOS_user_list3_upd(:,1),AOS_user_list3_upd(:,2),color,'filled','HandleVisibility','off');
    scatter(AOS_user_list3_upd(:,3),AOS_user_list3_upd(:,4),color,'filled','HandleVisibility','off');
    end
%     legend('CT','init DT','DT AOS','DT mAOS');
end

function [CUEs,CUE1_x, CUE1_y,CUE2_x, CUE2_y,CUE3_x, CUE3_y,...
    CT_BS_gain1,CT_BS_gain2,CT_BS_gain3]=generate_3_CTs(eNB1_x,eNB1_y,eNB2_x,eNB2_y,eNB3_x,eNB3_y,Cell_Radius)
        [CUE1_x, CUE1_y] = generate_CT(eNB1_x,eNB1_y,Cell_Radius);
        [CUE2_x, CUE2_y] = generate_CT(eNB2_x,eNB2_y,Cell_Radius);
        [CUE3_x, CUE3_y] = generate_CT(eNB3_x,eNB3_y,Cell_Radius);
        CT_BS_gain1 = LTE_channel_model_urban_micro_NLOS(CUE1_x,CUE1_y,eNB1_x,eNB1_y); 
        CT_BS_gain2 = LTE_channel_model_urban_micro_NLOS(CUE2_x,CUE2_y,eNB2_x,eNB2_y);
        CT_BS_gain3 = LTE_channel_model_urban_micro_NLOS(CUE3_x,CUE3_y,eNB3_x,eNB3_y);
        CUEs=[CUE1_x, CUE1_y;CUE2_x, CUE2_y;CUE3_x, CUE3_y];
    end 
function [AOS_user_list1_upd,AOS_user_list2_upd,AOS_user_list3_upd,...
    D2Ds_left1,D2Ds_left2,D2Ds_left3,...
    AOS_user_list1,AOS_user_list2,AOS_user_list3,...
    SE_new_AOS,SINR_C_new_all_AOS_dB,SINR_D_i_new_all_AOS_dB,numbers_of_pairs_new_AOS,...
    rank_AOS1,rank_AOS2,rank_AOS3]=...
    multi_cell_AOS...
    (user_list1,user_list2,user_list3,...
    Cell_Radius,D2D_Sep_Max,...
    Max_Users1,Max_Users2,Max_Users3,...
    CUE_Exp,DUE_Exp,CUE_SINR_min,DUE_SINR_min,...
    eNB1_x,eNB2_x,eNB3_x,eNB1_y,eNB2_y,eNB3_y,...
    CUE1_x, CUE1_y,CUE2_x, CUE2_y,CUE3_x, CUE3_y,...
    CT_BS_gain1,CT_BS_gain2,CT_BS_gain3) 
    
    SINR_D_i_new_all_AOS=[];
    SINR_C_new_all_AOS=[];

    [~,N_selected_PRS1,rank_AOS1,N_selected_AOS1]=single_cell_PRS_AOS(user_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,Max_Users1,Cell_Radius,D2D_Sep_Max);
%     disp('# of selected pairs (PRS/AOS/PC) in cell 1 = ');
    to_disp=[num2str(N_selected_PRS1),'/',num2str(N_selected_AOS1)];
%     disp(to_disp);
    [~,N_selected_PRS2,rank_AOS2,N_selected_AOS2]=single_cell_PRS_AOS(user_list2,eNB2_x,eNB2_y,CUE2_x,CUE2_y,Max_Users2,Cell_Radius,D2D_Sep_Max); 
%     disp('# of selected pairs (PRS/AOS/PC) in cell 2 = ');
    to_disp=[num2str(N_selected_PRS2),'/',num2str(N_selected_AOS2)];
%     disp(to_disp);
    [~,N_selected_PRS3,rank_AOS3,N_selected_AOS3]=single_cell_PRS_AOS(user_list3,eNB3_x,eNB3_y,CUE3_x,CUE3_y,Max_Users3,Cell_Radius,D2D_Sep_Max); 
%     disp('# of selected pairs (PRS/AOS/PC) in cell 3 = ');
    to_disp=[num2str(N_selected_PRS3),'/',num2str(N_selected_AOS3)];
%     disp(to_disp); 

    AOS_user_list1 = mask (user_list1,rank_AOS1); 
    AOS_user_list2 = mask (user_list2,rank_AOS2);
    AOS_user_list3 = mask (user_list3,rank_AOS3);

    [min_index12] = closest_pair_to_CT (AOS_user_list1,CUE2_x,CUE2_y);
    [min_index13] = closest_pair_to_CT (AOS_user_list1,CUE3_x,CUE3_y);
    [min_index21] = closest_pair_to_CT (AOS_user_list2,CUE1_x,CUE1_y);
    [min_index23] = closest_pair_to_CT (AOS_user_list2,CUE3_x,CUE3_y);
    [min_index31] = closest_pair_to_CT (AOS_user_list3,CUE1_x,CUE1_y);
    [min_index32] = closest_pair_to_CT (AOS_user_list3,CUE2_x,CUE2_y);

    %CT forbidden region calculation 
%     disp ' '; disp('CT check, cell 1 to 2:');
    [AOS_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (AOS_user_list1, min_index12,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,DUE_SINR_min,CUE_Exp,DUE_Exp,CUE3_x,CUE3_y); 
%     disp ' '; disp('CT check, cell 1 to 3:');                                                
    [AOS_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (AOS_user_list1_upd, min_index13,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,DUE_SINR_min,CUE_Exp,DUE_Exp,CUE2_x,CUE2_y); 
%     disp ' '; disp('CT check, cell 2 to 1:');                                                
    [AOS_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (AOS_user_list2, min_index21,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,DUE_SINR_min,CUE_Exp,DUE_Exp,CUE3_x,CUE3_y);
%     disp ' '; disp('CT check, cell 2 to 3:');                                                
    [AOS_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (AOS_user_list2_upd, min_index23,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,DUE_SINR_min,CUE_Exp,DUE_Exp,CUE1_x,CUE1_y); 
%     disp ' '; disp('CT check, cell 3 to 1:');                                                
    [AOS_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (AOS_user_list3, min_index31,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,DUE_SINR_min,CUE_Exp,DUE_Exp,CUE2_x,CUE2_y);
%     disp ' '; disp('CT check, cell 3 to 2:');                                                
    [AOS_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (AOS_user_list3_upd, min_index32,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,DUE_SINR_min,CUE_Exp,DUE_Exp,CUE1_x,CUE1_y); 
    %DT forbidden region calculation
%     disp ' '; disp('DT check, cell 1 to 2:'); %FR violated -> pair will be deleted from 1st list
    [AOS_user_list1_upd,min_index12,min_index13] = forbidden_region_DUE (AOS_user_list1_upd, min_index12, ...
                                                        AOS_user_list2_upd, min_index21,DUE_SINR_min,DUE_Exp, ...
                                                                        AOS_user_list3_upd); 
%     disp ' '; disp('DT check, cell 2 to 3:'); 
    [AOS_user_list2_upd,min_index23,min_index21] = forbidden_region_DUE (AOS_user_list2_upd, min_index23, ...
                                                        AOS_user_list3_upd, min_index32,DUE_SINR_min,DUE_Exp, ...
                                                                        AOS_user_list1_upd); 
%     disp ' '; disp('DT check, cell 3 to 1:'); 
    [AOS_user_list3_upd,min_index31,min_index32] = forbidden_region_DUE (AOS_user_list3_upd, min_index31, ...
                                                        AOS_user_list1_upd, min_index13,DUE_SINR_min,DUE_Exp, ...
                                                                        AOS_user_list2_upd); 
    n=1;
    numbers_of_pairs_new_AOS(n,1) = size(AOS_user_list1_upd,1)+size(AOS_user_list2_upd,1)+size(AOS_user_list3_upd,1);
        
    CT_user_list_all=[CUE1_x,CUE1_y;CUE2_x,CUE2_y;CUE3_x,CUE3_y];
    AOS_user_list_all_raw= [AOS_user_list1;AOS_user_list2;AOS_user_list3];
    AOS_user_list_all_new= [AOS_user_list1_upd;AOS_user_list2_upd;AOS_user_list3_upd];
        
    [SINR_C1_new_AOS,SINR_D_i1_new_AOS] = SINR_calc (1,AOS_user_list1_upd,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_new_AOS,SINR_D_i2_new_AOS] = SINR_calc (2,AOS_user_list2_upd,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_new_AOS,SINR_D_i3_new_AOS] = SINR_calc (3,AOS_user_list3_upd,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);

    mult=1;
    for k=1:size(SINR_D_i1_new_AOS,1)
        mult=mult * ( 1 + SINR_D_i1_new_AOS(k) );
    end
    for k=1:size(SINR_D_i2_new_AOS,1)
        mult=mult * ( 1 + SINR_D_i2_new_AOS(k) );
    end
    for k=1:size(SINR_D_i3_new_AOS,1)
        mult=mult * ( 1 + SINR_D_i3_new_AOS(k) );
    end
    SE_new_AOS(n)=log2( (1+SINR_C1_new_AOS)*(1+SINR_C2_new_AOS)*(1+SINR_C3_new_AOS) * mult );
    
    SINR_C_new_all_AOS=[SINR_C_new_all_AOS;SINR_C1_new_AOS;SINR_C2_new_AOS;SINR_C3_new_AOS];
    SINR_D_i_new_all_AOS=[SINR_D_i_new_all_AOS;SINR_D_i1_new_AOS;SINR_D_i2_new_AOS;SINR_D_i3_new_AOS];
    SINR_C_new_all_AOS_dB=pow2db(SINR_C_new_all_AOS); %multiAOS
    SINR_D_i_new_all_AOS_dB=pow2db(SINR_D_i_new_all_AOS);

    D2Ds_left1 = setdiff(user_list1,AOS_user_list1_upd,'rows');
    D2Ds_left2 = setdiff(user_list2,AOS_user_list2_upd,'rows');
    D2Ds_left3 = setdiff(user_list3,AOS_user_list3_upd,'rows');
end

function [ans] = dist (x1, y1, x2, y2) 
    ans= sqrt((abs(x1-x2))^2+(abs(y1-y2))^2);
end
function [SINR_C,SINR_D_i] = SINR_calc (numb_cell,user_list,eNB_x,eNB_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,user_list_all,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3)
    if (isempty(user_list))
        SINR_C=0;
        SINR_D_i=0;
        return
    end
    Bandwidth_kHz = 180; 
    Noise_dBm =	-174; 
    Noise_Total_Watts =	Bandwidth_kHz*10^3*(db2pow(Noise_dBm-30)); %for SINR calc
    CT_Tx_Power = 50*10^(-3); 
    %CT
    if numb_cell == 1 
        CT_BS_gain=CT_BS_gain1; %CT-BS gain in current cell, urban channel
        CT_BS_gain_other_1 = LTE_channel_model_urban_micro_NLOS(CUE2_x,CUE2_y,eNB_x,eNB_y); %other CT-current BS gain
        CT_BS_gain_other_2 = LTE_channel_model_urban_micro_NLOS(CUE3_x,CUE3_y,eNB_x,eNB_y); %other CT-current BS gain
    end
    if numb_cell == 2 
        CT_BS_gain=CT_BS_gain2;
        CT_BS_gain_other_1 = LTE_channel_model_urban_micro_NLOS(CUE1_x,CUE1_y,eNB_x,eNB_y); 
        CT_BS_gain_other_2 = LTE_channel_model_urban_micro_NLOS(CUE3_x,CUE3_y,eNB_x,eNB_y);
    end
    if numb_cell == 3 
        CT_BS_gain=CT_BS_gain3;
        CT_BS_gain_other_1 = LTE_channel_model_urban_micro_NLOS(CUE1_x,CUE1_y,eNB_x,eNB_y); 
        CT_BS_gain_other_2 = LTE_channel_model_urban_micro_NLOS(CUE2_x,CUE2_y,eNB_x,eNB_y);
    end
    
    %CT
    sum1=0; %sum over DTs' interference
    P_DT_i = 1*10^(-3);
    for i=1:size(user_list_all,1)
        DT_eNB_gain = LTE_channel_model_urban_micro_NLOS(user_list_all(i,1),user_list_all(i,2),eNB_x,eNB_y);%urban: channel model DT-eNB
        sum1=sum1+ (P_DT_i * DT_eNB_gain); %sum over DTs' interference towards BS
    end
    sum2=CT_Tx_Power*CT_BS_gain_other_1 + CT_Tx_Power*CT_BS_gain_other_2; %sum over other CTs' interference
    SINR_C = (CT_Tx_Power * CT_BS_gain) / (Noise_Total_Watts + sum1 ); %sum2 affect is minor
    
    %DT
    P_DT_i = 1*10^(-3); 
    for i=1:size(user_list,1)
        sum3=0; %sum over DTs' interference %%%was out of the for loop in ARP
        DT_Pair_gain_i = LTE_channel_model_indoor_hotspot_NLOS(user_list(i,1),user_list(i,2),user_list(i,3),user_list(i,4));
        G_CD_i1 = LTE_channel_model_indoor_hotspot_NLOS(user_list(i,3),user_list(i,4),CUE1_x,CUE1_y); %in Abu's code there is indoor                                                                              
        G_CD_i2 = LTE_channel_model_indoor_hotspot_NLOS(user_list(i,3),user_list(i,4),CUE2_x,CUE2_y);
        G_CD_i3 = LTE_channel_model_indoor_hotspot_NLOS(user_list(i,3),user_list(i,4),CUE3_x,CUE3_y);
        for j=1:size(user_list_all,1)
            if ( user_list(i,1) ~= user_list_all(j,1) )
                P_DT_j = 1*10^(-3); 
                G_D_ji = LTE_channel_model_indoor_hotspot_NLOS(user_list_all(j,1), ... %%%was urban in ARP, in Abu's code there is indoor
                    user_list_all(j,2),user_list(i,3),user_list(i,4));
                sum3=sum3+(P_DT_j*G_D_ji);
            end       
        end
        sum4=CT_Tx_Power*G_CD_i1+CT_Tx_Power*G_CD_i2+CT_Tx_Power*G_CD_i3; %sum over CTs' interference
        SINR_D_i(i) = (P_DT_i * DT_Pair_gain_i)/(sum3 + sum4 + Noise_Total_Watts);
    end           
    SINR_D_i=transpose(SINR_D_i);
end

%fct chooses rows from D2D_user_list with mask rank_PRS (w/ positive values)
function PRS_user_list = mask ( D2D_user_list, rank_PRS) 
    k=1;
    for i=1:size(rank_PRS,1)
       if rank_PRS(i)>0
          PRS_user_list(k,:) = D2D_user_list(i,:);
          k=k+1;
       end
    end
    if (k==1) 
        PRS_user_list=[];
    end
end

function [min_index] = closest_pair_to_CT (user_list, CUE_x, CUE_y)
   if (isempty(user_list))
       min_index=-1;
   else
       for i=1:size(user_list,1)
           pairs_list(i,1)=(user_list(i,1)+user_list(i,3))/2;
           pairs_list(i,2)=(user_list(i,2)+user_list(i,4))/2;
       end
       for i=1:size(pairs_list,1)
           dist_matrix(i)=dist(pairs_list(i,1),pairs_list(i,2),CUE_x,CUE_y); %(x1, y1, x2, y2)
       end
       [min_val,min_index]=min(dist_matrix);
   end
end

%fct calculates forbidden region CT/DT
function [PRS_user_list, min_index, min_index_to_other_cell] = forbidden_region_CUE (PRS_user_list, min_index, ...
                                        CUE_x, CUE_y,eNB_x,eNB_y,SINR_C,SINR_D,CUE_Exp,DUE_Exp,CUE_x_other,CUE_y_other) 
    while 1  
        if (isempty(PRS_user_list))  %inner cell list empty? -> neg indexes, fct ends
%             disp('list is empty');
            min_index=-1;
            min_index_to_other_cell=-1;
            break
        end 
        l=dist(PRS_user_list(min_index,1),PRS_user_list(min_index,2), ... 
            PRS_user_list(min_index,3),PRS_user_list(min_index,4));            %Tx DT- Rx DT
        r0=dist(CUE_x,CUE_y,eNB_x,eNB_y);                                     %CT-BS
        r1=dist(PRS_user_list(min_index,1),PRS_user_list(min_index,2), ... 
            eNB_x,eNB_y);                                                   %%% Tx DT - BS
        R=(SINR_C ^(1/DUE_Exp))*(SINR_D ^(1/CUE_Exp)) * r0^(CUE_Exp/DUE_Exp) * l / r1^(CUE_Exp/DUE_Exp);
        if R >= r1 %forbidden region is violated
            PRS_user_list(min_index,:)=[];
            [min_index] = closest_pair_to_CT (PRS_user_list,CUE_x,CUE_y); %after deleting pair, closest pairs in both directions should be recalculated
            [min_index_to_other_cell] = closest_pair_to_CT (PRS_user_list,CUE_x_other,CUE_y_other);
%             disp('pair deleted');
        else
            [min_index] = closest_pair_to_CT (PRS_user_list,CUE_x,CUE_y); 
            [min_index_to_other_cell] = closest_pair_to_CT (PRS_user_list,CUE_x_other,CUE_y_other);
            break
        end   
    end
end

function [min_index] = closest_pairs_d2d (user_list1, user_list2) 
   if ( (isempty(user_list1)) || (isempty(user_list2)) )
       min_index=-1;
   else
       for i=1:size(user_list1,1)
           pairs_list1(i,1)=(user_list1(i,1)+user_list1(i,3))/2;
           pairs_list1(i,2)=(user_list1(i,2)+user_list1(i,4))/2;
       end
       for i=1:size(user_list2,1)
           pairs_list2(i,1)=(user_list2(i,1)+user_list2(i,3))/2;
           pairs_list2(i,2)=(user_list2(i,2)+user_list2(i,4))/2;
       end
       for i=1:size(pairs_list1,1)
           for j=1:size(pairs_list2,1)
               dist_matrix(i,j)=dist(pairs_list1(i,1),pairs_list1(i,2),pairs_list2(j,1),pairs_list2(j,2)); %(x1, y1, x2, y2)
           end
       end
       minimum = min(min(dist_matrix)); %finding the min value in matrix, taking its x & y
       [row,col]=find(dist_matrix==minimum);%row-1st list,col-2nd list
       min_index=row;
   end
end

%fct calculates forbidden region DT/DT
function [PRS_user_list1, min_index1, min_index_to_other_cell] = forbidden_region_DUE (PRS_user_list1, min_index1, ...
                                        PRS_user_list2, min_index2,SINR_D,DUE_Exp,PRS_user_list_other) 
    while 1 
        if (isempty(PRS_user_list1)) %inner cell list empty? -> fct returns neg flags
%             disp('inner cell is D2D-pair-empty');
            min_index1=-1;
            min_index_to_other_cell=-1;
            break
        end 
        if min_index2==-1 %outer cell list empty? -> fct returns input values
%             disp('outer cell is D2D-pair-empty');
            % min_index1 stays the same
            [min_index_to_other_cell] = closest_pairs_d2d (PRS_user_list1,PRS_user_list_other);
            break
        end
        l=dist(PRS_user_list1(min_index1,1),PRS_user_list1(min_index1,2), ... %Tx DT- Rx DT of pair1:can be deleted
            PRS_user_list1(min_index1,3),PRS_user_list1(min_index1,4)); 
        r0=dist(PRS_user_list2(min_index2,1),PRS_user_list2(min_index2,2), ... %Tx DT- Rx DT of pair2
            PRS_user_list2(min_index2,3),PRS_user_list2(min_index2,4)); 
        r1=dist(PRS_user_list1(min_index1,1),PRS_user_list1(min_index1,2), ... %pair1 Tx - pair2 Rx 
            PRS_user_list2(min_index2,3),PRS_user_list2(min_index2,4)); 
        R= (SINR_D ^(2/DUE_Exp)) * l * r0 / r1; %R is around pair1
        if R >= r1 %forbidden region is violated
            PRS_user_list1(min_index1,:)=[];
            [min_index1] = closest_pairs_d2d (PRS_user_list1,PRS_user_list2); %after deleting pair, closest pairs in both directions should be recalculated
            [min_index_to_other_cell] = closest_pairs_d2d (PRS_user_list1,PRS_user_list_other);
%             disp('pair deleted');
        else
            [min_index1] = closest_pairs_d2d (PRS_user_list1,PRS_user_list2);
            [min_index_to_other_cell] = closest_pairs_d2d (PRS_user_list1,PRS_user_list_other);
            break
        end   
    end
end

function [D2D_user_list,N_selected_PC] = single_cell_PC_Abu_channel(D2D_user_list, eNB_x,eNB_y,CUE_x,CUE_y)
    SINR_th_BS=10; %dB
    SINR_th_DT=19; %dB
    SINR_th_BS_W=db2pow(SINR_th_BS); %6dB = 3.9811 Watts
    SINR_th_DT_W=db2pow(SINR_th_DT); %15dB = 31.6228 Watts

    %sensitivity of receivers is probably -100dBm, ignore it 
    
    Pd = 1*10^(-3); %W
    Pc = 50*10^(-3); %W
    Bandwidth_kHz = 180; 
    Noise_dBm =	-174; 
    Noise_Total_Watts =	Bandwidth_kHz*10^3*(db2pow(Noise_dBm-30)); 
    while 1
        %CT
        DT_rx_power=[];
        sum1=0; %sum over DTs interference towards BS      
        for i=1:size(D2D_user_list,1)
            sum1= sum1 + (Pd * LTE_channel_model_urban_micro_NLOS(D2D_user_list(i,1),D2D_user_list(i,2),eNB_x,eNB_y)); 
        end
        BS_rx_power = Pc * LTE_channel_model_urban_micro_NLOS(CUE_x,CUE_y,eNB_x,eNB_y) / (sum1 + Noise_Total_Watts);
        
        %DT
        for i=1:size(D2D_user_list,1)
            sum2=0; %sum over DTs interference towards i-th DTrx    
            for j=1:size(D2D_user_list,1)
                if (j~=i)
                    sum2= sum2 + Pd * LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(j,1),D2D_user_list(j,2),D2D_user_list(i,3),D2D_user_list(i,4));           
                end
            end
            DT_rx_power(i) = Pd * LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(i,1),D2D_user_list(i,2),D2D_user_list(i,3),D2D_user_list(i,4)) / ...
                        ( Pc * LTE_channel_model_indoor_hotspot_NLOS(CUE_x,CUE_y,D2D_user_list(i,3),D2D_user_list(i,4)) + sum2 + Noise_Total_Watts);  %was urban
        end 
        
        DT_rx_power=transpose(DT_rx_power);
        minimum = min(DT_rx_power); 
        row=find(DT_rx_power==minimum);
        if ( (isempty(D2D_user_list))  || ( (BS_rx_power>SINR_th_BS_W) &&  (DT_rx_power(row)>SINR_th_DT_W) ) ) %
            N_selected_PC=size(D2D_user_list,1);
            break
        else
            D2D_user_list(row,:)=[];
        end
    end
end

function [CUE_x,CUE_y]=generate_CT(eNB_x,eNB_y,Cell_Radius)       
    UE_Dist_Min = 10; % Minimum distance of any UE (i.e. CUE or DUE) from either the BS or another UE (i.e. CUE or DUE)
    locUE = UE_Dist_Min + (Cell_Radius - UE_Dist_Min)*sqrt(rand(1,1));
    theta_= 2*pi*rand(1,1);
    CUE_x = locUE*cos(theta_) + eNB_x ;
    CUE_y = locUE*sin(theta_) + eNB_y ;
end
