clc;clear;close all;format compact;

SINR_D_i_rand_all=[];
SINR_C_rand_all=[];
SINR_D_i_raw_all=[];
SINR_C_raw_all=[];
SINR_D_i_new_all=[];
SINR_C_new_all=[];
SINR_D_i_raw_all_AOS=[];
SINR_C_raw_all_AOS=[];
SINR_D_i_new_all_AOS=[];
SINR_C_new_all_AOS=[];
SINR_D_i_raw_all_PC=[];
SINR_C_raw_all_PC=[];
SINR_D_i_new_all_PC=[];
SINR_C_new_all_PC=[];

slice=1; %-not used- counter for different SINR thresholds
dB_value_grid=[33]; % CT SINR th, level of forb reg harshness
DUE_SINR_min_dB	= 46; %DT SINR th     6/15 ARP, 30/46 opt
for dB_value=[dB_value_grid] %-not used- loop over different SINR thresholds
    
no_runs=1000;
for n=1:no_runs
    to_disp=['iteration # ',num2str(n)];
    disp(to_disp);
    Cell_Radius = 200;
    D2D_Sep_Max = 0.1*Cell_Radius;
    Max_Users = 100; % pairs
    CUE_Exp = 3.67; % Cellular Pathloss Exponent - Alpha_zero
    DUE_Exp = 4.33; % Devices Pathloss Exponent - Alpha_d
    CUE_SINR_min_dB	= dB_value; 
    CUE_SINR_min = db2pow(CUE_SINR_min_dB); 
    DUE_SINR_min = db2pow(DUE_SINR_min_dB); 

    %coordinates of BSs, 3 cells
    eNB1_x = 200; eNB1_y = 200; 
    eNB2_x = 400; eNB2_y = 200+200*sqrt(3);
    eNB3_x = 600; eNB3_y = 200;

    %distibution of users in 3 cells
    D2D_user_list1 = LTE_UE_uniform_distribution_upd(eNB1_x,eNB1_y,Cell_Radius,D2D_Sep_Max, Max_Users); %Tx_x,Tx_y,Rx_x,Rx_y
    D2D_user_list2 = LTE_UE_uniform_distribution_upd(eNB2_x,eNB2_y,Cell_Radius,D2D_Sep_Max, Max_Users);
    D2D_user_list3 = LTE_UE_uniform_distribution_upd(eNB3_x,eNB3_y,Cell_Radius,D2D_Sep_Max, Max_Users);
    
    [CUE1_x, CUE1_y] = generate_CT(eNB1_x,eNB1_y,Cell_Radius);
    [CUE2_x, CUE2_y] = generate_CT(eNB2_x,eNB2_y,Cell_Radius);
    [CUE3_x, CUE3_y] = generate_CT(eNB3_x,eNB3_y,Cell_Radius);
    %choose D2D pairs within single cell 
    [rank_PRS1,N_selected_PRS1,rank_AOS1,N_selected_AOS1]=single_cell_PRS_AOS(D2D_user_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,Max_Users,Cell_Radius);
%     [PC_user_list1,N_selected_PC1]=single_cell_PC_Abu_channel(D2D_user_list1,eNB1_x,eNB1_y,CUE1_x, CUE1_y); %PC - power control, choose best pairs according to their SINR
    disp('# of selected pairs (PRS/AOS/PC) in cell 1 = ');
    to_disp=[num2str(N_selected_PRS1),'/',num2str(N_selected_AOS1)];
    disp(to_disp);
    [rank_PRS2,N_selected_PRS2,rank_AOS2,N_selected_AOS2]=single_cell_PRS_AOS(D2D_user_list2,eNB2_x,eNB2_y,CUE2_x,CUE2_y,Max_Users,Cell_Radius); 
%     [PC_user_list2,N_selected_PC2]=single_cell_PC_Abu_channel(D2D_user_list2,eNB2_x,eNB2_y,CUE2_x, CUE2_y);
    disp('# of selected pairs (PRS/AOS/PC) in cell 2 = ');
    to_disp=[num2str(N_selected_PRS2),'/',num2str(N_selected_AOS2)];
    disp(to_disp);
    [rank_PRS3,N_selected_PRS3,rank_AOS3,N_selected_AOS3]=single_cell_PRS_AOS(D2D_user_list3,eNB3_x,eNB3_y,CUE3_x,CUE3_y,Max_Users,Cell_Radius); 
%     [PC_user_list3,N_selected_PC3]=single_cell_PC_Abu_channel(D2D_user_list3,eNB3_x,eNB3_y,CUE3_x, CUE3_y);
    disp('# of selected pairs (PRS/AOS/PC) in cell 3 = ');
    to_disp=[num2str(N_selected_PRS3),'/',num2str(N_selected_AOS3)];
    disp(to_disp);

    PRS_user_list1 = mask (D2D_user_list1,rank_PRS1); % PRS/AOS/PC uses its own SINR CT & DT thresholds (6&15)
    PRS_user_list2 = mask (D2D_user_list2,rank_PRS2);
    PRS_user_list3 = mask (D2D_user_list3,rank_PRS3);

    AOS_user_list1 = mask (D2D_user_list1,rank_AOS1); 
    AOS_user_list2 = mask (D2D_user_list2,rank_AOS2);
    AOS_user_list3 = mask (D2D_user_list3,rank_AOS3);

    %% -------------PRS-----------------------
    %calculation of closest pair to the border
    [min_index12] = closest_pair_to_CT (PRS_user_list1,CUE2_x,CUE2_y);
    [min_index13] = closest_pair_to_CT (PRS_user_list1,CUE3_x,CUE3_y);
    [min_index21] = closest_pair_to_CT (PRS_user_list2,CUE1_x,CUE1_y);
    [min_index23] = closest_pair_to_CT (PRS_user_list2,CUE3_x,CUE3_y);
    [min_index31] = closest_pair_to_CT (PRS_user_list3,CUE1_x,CUE1_y);
    [min_index32] = closest_pair_to_CT (PRS_user_list3,CUE2_x,CUE2_y);

    %CT forbidden region calculation 
    disp('--------multi-cell PRS--------');
    disp ' '; disp('CT check, cell 1 to 2:');
    [PRS_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (PRS_user_list1, min_index12,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y); 
    disp ' '; disp('CT check, cell 1 to 3:');                                                
    [PRS_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (PRS_user_list1_upd, min_index13,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y); 
    disp ' '; disp('CT check, cell 2 to 1:');                                                
    [PRS_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (PRS_user_list2, min_index21,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y);
    disp ' '; disp('CT check, cell 2 to 3:');                                                
    [PRS_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (PRS_user_list2_upd, min_index23,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
    disp ' '; disp('CT check, cell 3 to 1:');                                                
    [PRS_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (PRS_user_list3, min_index31,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y);
    disp ' '; disp('CT check, cell 3 to 2:');                                                
    [PRS_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (PRS_user_list3_upd, min_index32,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
    %DT forbidden region calculation
    disp ' '; disp('DT check, cell 1 to 2:'); %FR violated -> pair will be deleted from 1st list
    [PRS_user_list1_upd,min_index12,min_index13] = forbidden_region_DUE (PRS_user_list1_upd, min_index12, ...
                                                        PRS_user_list2_upd, min_index21,DUE_SINR_min,DUE_Exp, ...
                                                                        PRS_user_list3_upd); 
    disp ' '; disp('DT check, cell 2 to 3:'); 
    [PRS_user_list2_upd,min_index23,min_index21] = forbidden_region_DUE (PRS_user_list2_upd, min_index23, ...
                                                        PRS_user_list3_upd, min_index32,DUE_SINR_min,DUE_Exp, ...
                                                                        PRS_user_list1_upd); 
    disp ' '; disp('DT check, cell 3 to 1:'); 
    [PRS_user_list3_upd,min_index31,min_index32] = forbidden_region_DUE (PRS_user_list3_upd, min_index31, ...
                                                        PRS_user_list1_upd, min_index13,DUE_SINR_min,DUE_Exp, ...
                                                                        PRS_user_list2_upd); 
    %% -------------AOS-----------------------
    %calculation of closest pair to the border
    [min_index12] = closest_pair_to_CT (AOS_user_list1,CUE2_x,CUE2_y);
    [min_index13] = closest_pair_to_CT (AOS_user_list1,CUE3_x,CUE3_y);
    [min_index21] = closest_pair_to_CT (AOS_user_list2,CUE1_x,CUE1_y);
    [min_index23] = closest_pair_to_CT (AOS_user_list2,CUE3_x,CUE3_y);
    [min_index31] = closest_pair_to_CT (AOS_user_list3,CUE1_x,CUE1_y);
    [min_index32] = closest_pair_to_CT (AOS_user_list3,CUE2_x,CUE2_y);

    %CT forbidden region calculation 
    disp('--------multi-cell AOS--------');
    disp ' '; disp('CT check, cell 1 to 2:');
    [AOS_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (AOS_user_list1, min_index12,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y); 
    disp ' '; disp('CT check, cell 1 to 3:');                                                
    [AOS_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (AOS_user_list1_upd, min_index13,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y); 
    disp ' '; disp('CT check, cell 2 to 1:');                                                
    [AOS_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (AOS_user_list2, min_index21,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y);
    disp ' '; disp('CT check, cell 2 to 3:');                                                
    [AOS_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (AOS_user_list2_upd, min_index23,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
    disp ' '; disp('CT check, cell 3 to 1:');                                                
    [AOS_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (AOS_user_list3, min_index31,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y);
    disp ' '; disp('CT check, cell 3 to 2:');                                                
    [AOS_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (AOS_user_list3_upd, min_index32,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
    %DT forbidden region calculation
    disp ' '; disp('DT check, cell 1 to 2:'); %FR violated -> pair will be deleted from 1st list
    [AOS_user_list1_upd,min_index12,min_index13] = forbidden_region_DUE (AOS_user_list1_upd, min_index12, ...
                                                        AOS_user_list2_upd, min_index21,DUE_SINR_min,DUE_Exp, ...
                                                                        AOS_user_list3_upd); 
    disp ' '; disp('DT check, cell 2 to 3:'); 
    [AOS_user_list2_upd,min_index23,min_index21] = forbidden_region_DUE (AOS_user_list2_upd, min_index23, ...
                                                        AOS_user_list3_upd, min_index32,DUE_SINR_min,DUE_Exp, ...
                                                                        AOS_user_list1_upd); 
    disp ' '; disp('DT check, cell 3 to 1:'); 
    [AOS_user_list3_upd,min_index31,min_index32] = forbidden_region_DUE (AOS_user_list3_upd, min_index31, ...
                                                        AOS_user_list1_upd, min_index13,DUE_SINR_min,DUE_Exp, ...
                                                                        AOS_user_list2_upd); 
%     %% -------------PC-----------------------
%     %calculation of closest pair to the border
%     [min_index12] = closest_pair_to_CT (PC_user_list1,CUE2_x,CUE2_y);
%     [min_index13] = closest_pair_to_CT (PC_user_list1,CUE3_x,CUE3_y);
%     [min_index21] = closest_pair_to_CT (PC_user_list2,CUE1_x,CUE1_y);
%     [min_index23] = closest_pair_to_CT (PC_user_list2,CUE3_x,CUE3_y);
%     [min_index31] = closest_pair_to_CT (PC_user_list3,CUE1_x,CUE1_y);
%     [min_index32] = closest_pair_to_CT (PC_user_list3,CUE2_x,CUE2_y);
% 
%     %CT forbidden region calculation 
%     disp('--------PC--------');
%     disp ' '; disp('CT check, cell 1 to 2:');
%     [PC_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (PC_user_list1, min_index12,CUE2_x, ...
%                                                         CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y); 
%     disp ' '; disp('CT check, cell 1 to 3:');                                                
%     [PC_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (PC_user_list1_upd, min_index13,CUE3_x, ...
%                                                         CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y); 
%     disp ' '; disp('CT check, cell 2 to 1:');                                                
%     [PC_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (PC_user_list2, min_index21,CUE1_x, ...
%                                                         CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y);
%     disp ' '; disp('CT check, cell 2 to 3:');                                                
%     [PC_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (PC_user_list2_upd, min_index23,CUE3_x, ...
%                                                         CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
%     disp ' '; disp('CT check, cell 3 to 1:');                                                
%     [PC_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (PC_user_list3, min_index31,CUE1_x, ...
%                                                         CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y);
%     disp ' '; disp('CT check, cell 3 to 2:');                                                
%     [PC_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (PC_user_list3_upd, min_index32,CUE2_x, ...
%                                                         CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
%     %DT forbidden region calculation
%     disp ' '; disp('DT check, cell 1 to 2:'); %FR violated -> pair will be deleted from 1st list
%     [PC_user_list1_upd,min_index12,min_index13] = forbidden_region_DUE (PC_user_list1_upd, min_index12, ...
%                                                         PC_user_list2_upd, min_index21,DUE_SINR_min,DUE_Exp, ...
%                                                                         PC_user_list3_upd); 
%     disp ' '; disp('DT check, cell 2 to 3:'); 
%     [PC_user_list2_upd,min_index23,min_index21] = forbidden_region_DUE (PC_user_list2_upd, min_index23, ...
%                                                         PC_user_list3_upd, min_index32,DUE_SINR_min,DUE_Exp, ...
%                                                                         PC_user_list1_upd); 
%     disp ' '; disp('DT check, cell 3 to 1:'); 
%     [PC_user_list3_upd,min_index31,min_index32] = forbidden_region_DUE (PC_user_list3_upd, min_index31, ...
%                                                         PC_user_list1_upd, min_index13,DUE_SINR_min,DUE_Exp, ...
%                                                                         PC_user_list2_upd); 
    %% ---------------------------------------  
    random_list1=D2D_user_list1(1:15:end,:);
    random_list2=D2D_user_list2(1:15:end,:);
    random_list3=D2D_user_list3(1:15:end,:);
    random_list=[random_list1;random_list2;random_list3];
    numbers_of_pairs_rand(n,slice)=length(random_list1)+length(random_list2)+length(random_list3); 
    
    numbers_of_pairs_raw(n,slice) = length(PRS_user_list1)+length(PRS_user_list2)+length(PRS_user_list3); %before   
%     av_numbers_of_pairs_raw=mean(numbers_of_pairs_raw);
    numbers_of_pairs_new(n,slice) = length(PRS_user_list1_upd)+length(PRS_user_list2_upd)+length(PRS_user_list3_upd); %after
%     av_numbers_of_pairs_new=mean(numbers_of_pairs_new);
%     deleted_pairs (n,slice) = numbers_of_pairs_raw(n,slice) - numbers_of_pairs_new(n,slice); %difference 
%     av_deleted_pairs = mean(deleted_pairs);

    numbers_of_pairs_raw_AOS(n,slice) = length(AOS_user_list1)+length(AOS_user_list2)+length(AOS_user_list3); 
    numbers_of_pairs_new_AOS(n,slice) = length(AOS_user_list1_upd)+length(AOS_user_list2_upd)+length(AOS_user_list3_upd);
    
%     numbers_of_pairs_raw_PC(n,slice) = length(PC_user_list1)+length(PC_user_list2)+length(PC_user_list3); 
%     numbers_of_pairs_new_PC(n,slice) = length(PC_user_list1_upd)+length(PC_user_list2_upd)+length(PC_user_list3_upd);
    
    CT_user_list_all=[CUE1_x,CUE1_y;CUE2_x,CUE2_y;CUE3_x,CUE3_y];
    % D2D_user_list_all=[D2D_user_list1;D2D_user_list2;D2D_user_list3];
    PRS_user_list_all_raw= [PRS_user_list1;PRS_user_list2;PRS_user_list3];
    PRS_user_list_all_new= [PRS_user_list1_upd;PRS_user_list2_upd;PRS_user_list3_upd];
    AOS_user_list_all_raw= [AOS_user_list1;AOS_user_list2;AOS_user_list3];
    AOS_user_list_all_new= [AOS_user_list1_upd;AOS_user_list2_upd;AOS_user_list3_upd];
%     PC_user_list_all_raw= [PC_user_list1;PC_user_list2;PC_user_list3];
%     PC_user_list_all_new= [PC_user_list1_upd;PC_user_list2_upd;PC_user_list3_upd];
    
    CT_BS_gain1 = LTE_channel_model_urban_micro_NLOS(CUE1_x,CUE1_y,eNB1_x,eNB1_y); 
    CT_BS_gain2 = LTE_channel_model_urban_micro_NLOS(CUE2_x,CUE2_y,eNB2_x,eNB2_y);
    CT_BS_gain3 = LTE_channel_model_urban_micro_NLOS(CUE3_x,CUE3_y,eNB3_x,eNB3_y);

    %% -----------------SINR, SE random--------------------
    [SINR_C1_rand,SINR_D_i1_rand] = SINR_calc (1,random_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,random_list,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_rand,SINR_D_i2_rand] = SINR_calc (2,random_list2,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,random_list,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_rand,SINR_D_i3_rand] = SINR_calc (3,random_list3,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,random_list,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    mult=1;
    for k=1:length(SINR_D_i1_rand)
        mult=mult * ( 1 + SINR_D_i1_rand(k) );
    end
    for k=1:length(SINR_D_i2_rand)
        mult=mult * ( 1 + SINR_D_i2_rand(k) );
    end
    for k=1:length(SINR_D_i3_rand)
        mult=mult * ( 1 + SINR_D_i3_rand(k) );
    end
    SE_rand(n)=log2( (1+SINR_C1_rand)*(1+SINR_C2_rand)*(1+SINR_C3_rand) * mult );

    %% -----------------SINR, SE raw PRS--------------------
    [SINR_C1_raw,SINR_D_i1_raw] = SINR_calc (1,PRS_user_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PRS_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_raw,SINR_D_i2_raw] = SINR_calc (2,PRS_user_list2,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PRS_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_raw,SINR_D_i3_raw] = SINR_calc (3,PRS_user_list3,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PRS_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    mult=1;
    for k=1:length(SINR_D_i1_raw)
        mult=mult * ( 1 + SINR_D_i1_raw(k) );
    end
    for k=1:length(SINR_D_i2_raw)
        mult=mult * ( 1 + SINR_D_i2_raw(k) );
    end
    for k=1:length(SINR_D_i3_raw)
        mult=mult * ( 1 + SINR_D_i3_raw(k) );
    end
    SE_raw(n)=log2( (1+SINR_C1_raw)*(1+SINR_C2_raw)*(1+SINR_C3_raw) * mult );

    %% -----------------SINR, SE multicell PRS--------------------
    [SINR_C1_new,SINR_D_i1_new] = SINR_calc (1,PRS_user_list1_upd,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PRS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_new,SINR_D_i2_new] = SINR_calc (2,PRS_user_list2_upd,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PRS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_new,SINR_D_i3_new] = SINR_calc (3,PRS_user_list3_upd,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PRS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);

    mult=1;
    for k=1:length(SINR_D_i1_new)
        mult=mult * ( 1 + SINR_D_i1_new(k) );
    end
    for k=1:length(SINR_D_i2_new)
        mult=mult * ( 1 + SINR_D_i2_new(k) );
    end
    for k=1:length(SINR_D_i3_new)
        mult=mult * ( 1 + SINR_D_i3_new(k) );
    end
    SE_new(n)=log2( (1+SINR_C1_new)*(1+SINR_C2_new)*(1+SINR_C3_new) * mult );

    %% -----------------SINR, SE raw AOS--------------------
    [SINR_C1_raw_AOS,SINR_D_i1_raw_AOS] = SINR_calc (1,AOS_user_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_raw_AOS,SINR_D_i2_raw_AOS] = SINR_calc (2,AOS_user_list2,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_raw_AOS,SINR_D_i3_raw_AOS] = SINR_calc (3,AOS_user_list3,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    mult=1;
    for k=1:length(SINR_D_i1_raw_AOS)
        mult=mult * ( 1 + SINR_D_i1_raw_AOS(k) );
    end
    for k=1:length(SINR_D_i2_raw_AOS)
        mult=mult * ( 1 + SINR_D_i2_raw_AOS(k) );
    end
    for k=1:length(SINR_D_i3_raw_AOS)
        mult=mult * ( 1 + SINR_D_i3_raw_AOS(k) );
    end
    SE_raw_AOS(n)=log2( (1+SINR_C1_raw_AOS)*(1+SINR_C2_raw_AOS)*(1+SINR_C3_raw_AOS) * mult );

    %% -----------------SINR, SE multicell AOS--------------------
    [SINR_C1_new_AOS,SINR_D_i1_new_AOS] = SINR_calc (1,AOS_user_list1_upd,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_new_AOS,SINR_D_i2_new_AOS] = SINR_calc (2,AOS_user_list2_upd,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_new_AOS,SINR_D_i3_new_AOS] = SINR_calc (3,AOS_user_list3_upd,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,AOS_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);

    mult=1;
    for k=1:length(SINR_D_i1_new_AOS)
        mult=mult * ( 1 + SINR_D_i1_new_AOS(k) );
    end
    for k=1:length(SINR_D_i2_new_AOS)
        mult=mult * ( 1 + SINR_D_i2_new_AOS(k) );
    end
    for k=1:length(SINR_D_i3_new_AOS)
        mult=mult * ( 1 + SINR_D_i3_new_AOS(k) );
    end
    SE_new_AOS(n)=log2( (1+SINR_C1_new_AOS)*(1+SINR_C2_new_AOS)*(1+SINR_C3_new_AOS) * mult );
%     %% -----------------SINR, SE raw PC--------------------
%     [SINR_C1_raw_PC,SINR_D_i1_raw_PC] = SINR_calc (1,PC_user_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
%     [SINR_C2_raw_PC,SINR_D_i2_raw_PC] = SINR_calc (2,PC_user_list2,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
%     [SINR_C3_raw_PC,SINR_D_i3_raw_PC] = SINR_calc (3,PC_user_list3,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
% 
%     mult=1;
%     for k=1:length(SINR_D_i1_raw_PC)
%         mult=mult * ( 1 + SINR_D_i1_raw_PC(k) );
%     end
%     for k=1:length(SINR_D_i2_raw_PC)
%         mult=mult * ( 1 + SINR_D_i2_raw_PC(k) );
%     end
%     for k=1:length(SINR_D_i3_raw_PC)
%         mult=mult * ( 1 + SINR_D_i3_raw_PC(k) );
%     end
%     SE_raw_PC(n)=log2( (1+SINR_C1_raw_PC)*(1+SINR_C2_raw_PC)*(1+SINR_C3_raw_PC) * mult );
%     %% -----------------SINR, SE multicell PC--------------------
%     [SINR_C1_new_PC,SINR_D_i1_new_PC] = SINR_calc (1,PC_user_list1_upd,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
%     [SINR_C2_new_PC,SINR_D_i2_new_PC] = SINR_calc (2,PC_user_list1_upd,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
%     [SINR_C3_new_PC,SINR_D_i3_new_PC] = SINR_calc (3,PC_user_list1_upd,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
% 
%     mult=1;
%     for k=1:length(SINR_D_i1_new_PC)
%         mult=mult * ( 1 + SINR_D_i1_new_PC(k) );
%     end
%     for k=1:length(SINR_D_i2_new_PC)
%         mult=mult * ( 1 + SINR_D_i2_new_PC(k) );
%     end
%     for k=1:length(SINR_D_i3_new_PC)
%         mult=mult * ( 1 + SINR_D_i3_new_PC(k) );
%     end
%     SE_new_PC(n)=log2( (1+SINR_C1_new_PC)*(1+SINR_C2_new_PC)*(1+SINR_C3_new_PC) * mult );
    
    %% ---------------------------------------   
    SINR_C_rand_all=[SINR_C_rand_all;SINR_C1_rand;SINR_C2_rand;SINR_C3_rand];
    SINR_D_i_rand_all=[SINR_D_i_rand_all;SINR_D_i1_rand;SINR_D_i2_rand;SINR_D_i3_rand];
    
    SINR_C_raw_all=[SINR_C_raw_all;SINR_C1_raw;SINR_C2_raw;SINR_C3_raw];
    SINR_D_i_raw_all=[SINR_D_i_raw_all;SINR_D_i1_raw;SINR_D_i2_raw;SINR_D_i3_raw];
    SINR_C_new_all=[SINR_C_new_all;SINR_C1_new;SINR_C2_new;SINR_C3_new];
    SINR_D_i_new_all=[SINR_D_i_new_all;SINR_D_i1_new;SINR_D_i2_new;SINR_D_i3_new];
    
    SINR_C_raw_all_AOS=[SINR_C_raw_all_AOS;SINR_C1_raw_AOS;SINR_C2_raw_AOS;SINR_C3_raw_AOS];
    SINR_D_i_raw_all_AOS=[SINR_D_i_raw_all_AOS;SINR_D_i1_raw_AOS;SINR_D_i2_raw_AOS;SINR_D_i3_raw_AOS];
    SINR_C_new_all_AOS=[SINR_C_new_all_AOS;SINR_C1_new_AOS;SINR_C2_new_AOS;SINR_C3_new_AOS];
    SINR_D_i_new_all_AOS=[SINR_D_i_new_all_AOS;SINR_D_i1_new_AOS;SINR_D_i2_new_AOS;SINR_D_i3_new_AOS];
    
%     SINR_C_raw_all_PC=[SINR_C_raw_all_PC;SINR_C1_raw_PC;SINR_C2_raw_PC;SINR_C3_raw_PC];
%     SINR_D_i_raw_all_PC=[SINR_D_i_raw_all_PC;SINR_D_i1_raw_PC;SINR_D_i2_raw_PC;SINR_D_i3_raw_PC];
%     SINR_C_new_all_PC=[SINR_C_new_all_PC;SINR_C1_new_PC;SINR_C2_new_PC;SINR_C3_new_PC];
%     SINR_D_i_new_all_PC=[SINR_D_i_new_all_PC;SINR_D_i1_new_PC;SINR_D_i2_new_PC;SINR_D_i3_new_PC];
end   
% SINR_C_all_all_dB=pow2db(SINR_C_all_all);
% SINR_D_i_all_all_dB=pow2db(SINR_D_i_all_all);
SINR_C_rand_all_dB=pow2db(SINR_C_rand_all); %rand
SINR_D_i_rand_all_dB=pow2db(SINR_D_i_rand_all);

SINR_C_raw_all_dB=pow2db(SINR_C_raw_all); %PRS
SINR_D_i_raw_all_dB=pow2db(SINR_D_i_raw_all);
SINR_C_new_all_dB=pow2db(SINR_C_new_all); %multiPRS
SINR_D_i_new_all_dB=pow2db(SINR_D_i_new_all);

SINR_C_raw_all_AOS_dB=pow2db(SINR_C_raw_all_AOS); %AOS
SINR_D_i_raw_all_AOS_dB=pow2db(SINR_D_i_raw_all_AOS);
SINR_C_new_all_AOS_dB=pow2db(SINR_C_new_all_AOS); %multiAOS
SINR_D_i_new_all_AOS_dB=pow2db(SINR_D_i_new_all_AOS);

% SINR_C_raw_all_PC_dB=pow2db(SINR_C_raw_all_PC); %PC
% SINR_D_i_raw_all_PC_dB=pow2db(SINR_D_i_raw_all_PC);
% SINR_C_new_all_PC_dB=pow2db(SINR_C_new_all_PC); %multiPC
% SINR_D_i_new_all_PC_dB=pow2db(SINR_D_i_new_all_PC);

SE_rand=transpose(SE_rand);
SE_raw=transpose(SE_raw);
SE_new=transpose(SE_new);
SE_raw_AOS=transpose(SE_raw_AOS);
SE_new_AOS=transpose(SE_new_AOS);
% SE_raw_PC=transpose(SE_raw_PC);
% SE_new_PC=transpose(SE_new_PC);

% av_SE_rand(slice)=mean(SE_rand);
% av_SE_raw(slice)=mean(SE_raw);
% av_SE_new(slice)=mean(SE_new);
slice=slice+1;
end

%% PLOTS

% UT distribution visualization
figure
viscircles([200 200],200,'Color','k','linewidth',0.5);
viscircles([600 200],200,'Color','k','linewidth',0.5);
viscircles([400 200+200*sqrt(3)],200,'Color','k','linewidth',0.5);
axis ([-10 810 -10 800]);
axis square
grid on;
hold on;
scatter([CUE1_x CUE2_x CUE3_x],[CUE1_y CUE2_y CUE3_y],'green','filled');

scatter(D2D_user_list1(:,1),D2D_user_list1(:,2),'k','MarkerEdgeAlpha',.2);
scatter(D2D_user_list1(:,3),D2D_user_list1(:,4),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
scatter(D2D_user_list2(:,1),D2D_user_list2(:,2),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
scatter(D2D_user_list2(:,3),D2D_user_list2(:,4),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
scatter(D2D_user_list3(:,1),D2D_user_list3(:,2),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');
scatter(D2D_user_list3(:,3),D2D_user_list3(:,4),'k','MarkerEdgeAlpha',.2,'HandleVisibility','off');

% scatter(PC_user_list1(:,1),PC_user_list1(:,2),'blue');
% scatter(PC_user_list1(:,3),PC_user_list1(:,4),'blue','HandleVisibility','off');
% scatter(PC_user_list2(:,1),PC_user_list2(:,2),'blue','HandleVisibility','off');
% scatter(PC_user_list2(:,3),PC_user_list2(:,4),'blue','HandleVisibility','off');
% scatter(PC_user_list3(:,1),PC_user_list3(:,2),'blue','HandleVisibility','off');
% scatter(PC_user_list3(:,3),PC_user_list3(:,4),'blue','HandleVisibility','off');

scatter(AOS_user_list1(:,1),AOS_user_list1(:,2),'red');
scatter(AOS_user_list1(:,3),AOS_user_list1(:,4),'red','HandleVisibility','off');
scatter(AOS_user_list2(:,1),AOS_user_list2(:,2),'red','HandleVisibility','off');
scatter(AOS_user_list2(:,3),AOS_user_list2(:,4),'red','HandleVisibility','off');
scatter(AOS_user_list3(:,1),AOS_user_list3(:,2),'red','HandleVisibility','off');
scatter(AOS_user_list3(:,3),AOS_user_list3(:,4),'red','HandleVisibility','off');

% mixed1=intersect(PC_user_list1,AOS_user_list1,'rows'); %selected by both PC&AOS
% mixed2=intersect(PC_user_list2,AOS_user_list2,'rows');
% mixed3=intersect(PC_user_list3,AOS_user_list3,'rows');
% scatter(mixed1(:,1),mixed1(:,2),'magenta');
% scatter(mixed1(:,3),mixed1(:,4),'magenta','HandleVisibility','off');
% scatter(mixed2(:,1),mixed2(:,2),'magenta','HandleVisibility','off');
% scatter(mixed2(:,3),mixed2(:,4),'magenta','HandleVisibility','off');
% scatter(mixed3(:,1),mixed3(:,2),'magenta','HandleVisibility','off');
% scatter(mixed3(:,3),mixed3(:,4),'magenta','HandleVisibility','off');

% scatter(PC_user_list1_upd(:,1),PC_user_list1_upd(:,2),'blue','filled');
% scatter(PC_user_list1_upd(:,3),PC_user_list1_upd(:,4),'blue','filled','HandleVisibility','off');
% scatter(PC_user_list2_upd(:,1),PC_user_list2_upd(:,2),'blue','filled','HandleVisibility','off'); %HandleVisibility is off for legend to have only 3 opaque properties
% scatter(PC_user_list2_upd(:,3),PC_user_list2_upd(:,4),'blue','filled','HandleVisibility','off');
% scatter(PC_user_list3_upd(:,1),PC_user_list3_upd(:,2),'blue','filled','HandleVisibility','off');
% % scatter(PC_user_list3_upd(:,3),PC_user_list3_upd(:,4),'blue','filled','HandleVisibility','off');

scatter(AOS_user_list1_upd(:,1),AOS_user_list1_upd(:,2),'red','filled');
scatter(AOS_user_list1_upd(:,3),AOS_user_list1_upd(:,4),'red','filled','HandleVisibility','off');
scatter(AOS_user_list2_upd(:,1),AOS_user_list2_upd(:,2),'red','filled','HandleVisibility','off'); %HandleVisibility is off for legend to have only 3 opaque properties
scatter(AOS_user_list2_upd(:,3),AOS_user_list2_upd(:,4),'red','filled','HandleVisibility','off');
scatter(AOS_user_list3_upd(:,1),AOS_user_list3_upd(:,2),'red','filled','HandleVisibility','off');
scatter(AOS_user_list3_upd(:,3),AOS_user_list3_upd(:,4),'red','filled','HandleVisibility','off');

% mixed1_filled=intersect(PC_user_list1_upd,AOS_user_list1_upd,'rows');
% mixed2_filled=intersect(PC_user_list2_upd,AOS_user_list2_upd,'rows');
% mixed3_filled=intersect(PC_user_list3_upd,AOS_user_list3_upd,'rows');
% scatter(mixed1_filled(:,1),mixed1_filled(:,2),'magenta','filled');
% scatter(mixed1_filled(:,3),mixed1_filled(:,4),'magenta','filled','HandleVisibility','off');
% scatter(mixed2_filled(:,1),mixed2_filled(:,2),'magenta','filled','HandleVisibility','off');
% scatter(mixed2_filled(:,3),mixed2_filled(:,4),'magenta','filled','HandleVisibility','off');
% scatter(mixed3_filled(:,1),mixed3_filled(:,2),'magenta','filled','HandleVisibility','off');
% scatter(mixed3_filled(:,3),mixed3_filled(:,4),'magenta','filled','HandleVisibility','off');

% legend('CT','init DT','DT PC','DT AOS','DT PC&AOS','DT multiPC','DT multiAOS','DT multiPC&AOS');
legend('CT','init DT','DT AOS','DT multiAOS');


%CDF DT
figure
cdfplot(SINR_D_i_rand_all_dB);
xlim([0 50]);
hold on
% cdfplot(SINR_D_i_raw_all_dB);
cdfplot(SINR_D_i_new_all_dB);
% cdfplot(SINR_D_i_raw_all_AOS_dB);
cdfplot(SINR_D_i_new_all_AOS_dB);
% cdfplot(SINR_D_i_raw_all_PC_dB);
% cdfplot(SINR_D_i_new_all_PC_dB);
grid on
% legend('rand','PRS','multiPRS','AOS','multiAOS','PC','multiPC');
legend('random','multiPRS','multiAOS');
xlabel('DT SINR (dB)','FontName','Arial','FontSize',14);
ylabel('CDF','FontName','Arial','FontSize',14);

%CDF CT
figure
cdfplot(SINR_C_rand_all_dB);
xlim([0 50]);
hold on
% cdfplot(SINR_C_raw_all_dB);
cdfplot(SINR_C_new_all_dB);
% cdfplot(SINR_C_raw_all_AOS_dB);
cdfplot(SINR_C_new_all_AOS_dB);
% cdfplot(SINR_C_raw_all_PC_dB);
% cdfplot(SINR_C_new_all_PC_dB);
grid on
% legend('rand','PRS','multiPRS','AOS','multiAOS','PC','multiPC');
legend('random','multiPRS','multiAOS');
xlabel('CT SINR (dB)','FontName','Arial','FontSize',14);
ylabel('CDF','FontName','Arial','FontSize',14);

%SE
A=[numbers_of_pairs_rand,SE_rand]; 
B=[numbers_of_pairs_raw,SE_raw];
C=[numbers_of_pairs_new,SE_new];
D=[numbers_of_pairs_raw_AOS,SE_raw_AOS];
E=[numbers_of_pairs_new_AOS,SE_new_AOS];
% F=[numbers_of_pairs_raw_PC,SE_raw_PC];
% G=[numbers_of_pairs_new_PC,SE_new_PC];
A=sortrows(A); %from smallest # of pairs to highest
[Alink,aa] = findgroups(A(:,1)); %aa - different values of pairs, Alink - links from aa to A
A1 = [ aa, splitapply(@mean,A(:,2),Alink)]; %find mean value for each # of pairs
B=sortrows(B);
[Blink,bb] = findgroups(B(:,1));
B1 = [ bb, splitapply(@mean,B(:,2),Blink)];
C= sortrows(C);
[Clink,cc] = findgroups(C(:,1));
C1 = [ cc, splitapply(@mean,C(:,2),Clink)];
D= sortrows(D);
[Dlink,dd] = findgroups(D(:,1));
D1 = [ dd, splitapply(@mean,D(:,2),Dlink)];
E= sortrows(E);
[Elink,ee] = findgroups(E(:,1));
E1 = [ ee, splitapply(@mean,E(:,2),Elink)];
% F= sortrows(F);
% [Flink,ff] = findgroups(F(:,1));
% F1 = [ ff, splitapply(@mean,F(:,2),Flink)];
% G= sortrows(G);
% [Glink,gg] = findgroups(G(:,1));
% G1 = [ gg, splitapply(@mean,G(:,2),Glink)];

figure
plot(A1(:,1),A1(:,2),'b-o','linewidth',2.5); %b r y m g c 
hold on
% plot(B1(:,1),B1(:,2),'r-o','linewidth',2.5);
plot(C1(:,1),C1(:,2),'r-o','linewidth',2.5);
% plot(D1(:,1),D1(:,2),'m-o','linewidth',2.5);
plot(E1(:,1),E1(:,2),'y-o','linewidth',2.5);
% plot(F1(:,1),F1(:,2),'c-o','linewidth',2.5);
% plot(G1(:,1),G1(:,2),'m-o','linewidth',2.5);
grid on
% legend('rand','PRS','multiPRS','AOS','multiAOS','PC','multiPC');
legend('rand','multiPRS','multiAOS');
xlabel('Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);

save('A1.mat','A1');save('C1.mat','C1');save('E1.mat','E1');

%% FUNCTIONS
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
    for i=1:length(rank_PRS)
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
                                        CUE_x, CUE_y,eNB_x,eNB_y,SINR,Exp,CUE_x_other,CUE_y_other) 
    while 1  
        if (isempty(PRS_user_list))  %inner cell list empty? -> neg indexes, fct ends
            disp('list is empty');
            min_index=-1;
            min_index_to_other_cell=-1;
            break
        end 
        l=dist(PRS_user_list(min_index,1),PRS_user_list(min_index,2), ... 
            PRS_user_list(min_index,3),PRS_user_list(min_index,4));            %Tx DT- Rx DT
        r0=dist(CUE_x,CUE_y,eNB_x,eNB_y);                                     %CT-BS
        r1=dist(PRS_user_list(min_index,1),PRS_user_list(min_index,2), ... 
            eNB_x,eNB_y);                                                   %%% Tx DT - BS; was Rx DT-CT in ARP
        R= (SINR ^(2/Exp)) * l * r0 / r1;
        if R >= r1 %forbidden region is violated
            PRS_user_list(min_index,:)=[];
            [min_index] = closest_pair_to_CT (PRS_user_list,CUE_x,CUE_y); %after deleting pair, closest pairs in both directions should be recalculated
            [min_index_to_other_cell] = closest_pair_to_CT (PRS_user_list,CUE_x_other,CUE_y_other);
            disp('pair deleted');
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
                                        PRS_user_list2, min_index2,SINR,Exp,PRS_user_list_other) 
    while 1 
        if (isempty(PRS_user_list1)) %inner cell list empty? -> fct returns neg flags
            disp('inner cell is D2D-pair-empty');
            min_index1=-1;
            min_index_to_other_cell=-1;
            break
        end 
        if min_index2==-1 %outer cell list empty? -> fct returns input values
            disp('outer cell is D2D-pair-empty');
            % min_index1 stays the same
            [min_index_to_other_cell] = closest_pairs_d2d (PRS_user_list1,PRS_user_list_other);
            break
        end
        l=dist(PRS_user_list1(min_index1,1),PRS_user_list1(min_index1,2), ... %Tx DT- Rx DT of pair1:can be deleted
            PRS_user_list1(min_index1,3),PRS_user_list1(min_index1,4)); 
        r0=dist(PRS_user_list2(min_index2,1),PRS_user_list2(min_index2,2), ... %Tx DT- Rx DT of pair2
            PRS_user_list2(min_index2,3),PRS_user_list2(min_index2,4)); 
        r1=dist(PRS_user_list1(min_index1,1),PRS_user_list1(min_index1,2), ... %pair1 Tx - pair2 Rx ;was pair1 Rx - pair2 Rx in ARP
            PRS_user_list2(min_index2,3),PRS_user_list2(min_index2,4)); 
        R= (SINR ^(2/Exp)) * l * r0 / r1;
        if R >= r1 %forbidden region is violated
            PRS_user_list1(min_index1,:)=[];
            [min_index1] = closest_pairs_d2d (PRS_user_list1,PRS_user_list2); %after deleting pair, closest pairs in both directions should be recalculated
            [min_index_to_other_cell] = closest_pairs_d2d (PRS_user_list1,PRS_user_list_other);
            disp('pair deleted');
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
%% BACKUP
