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

dB_value=[33]; % CT SINR th, level of forb reg harshness
DUE_SINR_min_dB	= 46; %DT SINR th     6/15 ARP, 30/46 opt
grid_Kazan=xlsread('Points.xlsx');

no_runs=1;
for n=1:no_runs
    to_disp=['iteration # ',num2str(n)];
    disp(to_disp);
    Cell_Radius = 100;
    D2D_Sep_Max = 0.1*Cell_Radius;
    Max_Users = 50;
    CUE_Exp = 3.67; % Cellular Pathloss Exponent - Alpha_zero
    DUE_Exp = 4.33; % Devices Pathloss Exponent - Alpha_d
    CUE_SINR_min_dB	= dB_value; 
    CUE_SINR_min = db2pow(CUE_SINR_min_dB); 
    DUE_SINR_min = db2pow(DUE_SINR_min_dB); 
    eNB1_x = -272; eNB1_y = -645; 
    eNB2_x = -374; eNB2_y = -809;
    eNB3_x = -491; eNB3_y = -662;

    eNB_img_x=-384; eNB_img_y = -697; 
    Cell_Radius_img= 220;
    Max_Users_single = 100;
    
    load('C:\Users\Xiaomi\Documents\YASMP\ARP\codes\case3\D2D_user_list_all_mapped.mat','D2D_user_list_all_mapped');
    load('C:\Users\Xiaomi\Documents\YASMP\ARP\codes\case3\CT_user_list_all_mapped.mat','CT_user_list_all_mapped'); 
    CUE1_x=CT_user_list_all_mapped(1,1);    CUE1_y=CT_user_list_all_mapped(1,2);
    CUE2_x=CT_user_list_all_mapped(2,1);    CUE2_y=CT_user_list_all_mapped(2,2);
    CUE3_x=CT_user_list_all_mapped(3,1);    CUE3_y=CT_user_list_all_mapped(3,2);
    
    Loss_dB=readtable('C:\Users\Xiaomi\Documents\YASMP\ARP\codes\case3\case3.xls', 'Sheet', 3); 
    %1st row in .xls should be empty to be read properly
    Loss_dB=table2array(Loss_dB); 
    for i=1:size(Loss_dB,1)
        for j=1:size(Loss_dB,2)
            Loss(i,j)=10^(Loss_dB(i,j)/10); %times
            Gain(i,j)=1/Loss(i,j);
        end
    end
    
    %assigning initialized DTs to BSs
    [D2D_user_list1,D2D_user_list2,D2D_user_list3]=assign_to_cell_via_Gain(D2D_user_list_all_mapped,Gain);

    Max_Users1=length(D2D_user_list1); Max_Users2=length(D2D_user_list2); Max_Users3=length(D2D_user_list3);
    
    %choose D2D pairs within single cell 
    [rank_PRS1,N_selected_PRS1,rank_AOS1,N_selected_AOS1]=single_cell_PRS_AOS(D2D_user_list1,eNB1_x,eNB1_y,CUE1_x, CUE1_y,Max_Users1,Cell_Radius);
    [PC_user_list1,N_selected_PC1]=single_cell_PC_Abu_channel(D2D_user_list1,eNB1_x,eNB1_y,CUE1_x, CUE1_y); %PC - power control, choose best pairs according to their SINR
    disp('# of selected pairs (PRS/AOS/PC) in cell 1 = ');
    to_disp=[num2str(N_selected_PRS1),'/',num2str(N_selected_AOS1),'/',num2str(N_selected_PC1)];
    disp(to_disp);
    [rank_PRS2,N_selected_PRS2,rank_AOS2,N_selected_AOS2]=single_cell_PRS_AOS(D2D_user_list2,eNB2_x,eNB2_y,CUE2_x, CUE2_y,Max_Users2,Cell_Radius); 
    [PC_user_list2,N_selected_PC2]=single_cell_PC_Abu_channel(D2D_user_list2,eNB2_x,eNB2_y,CUE2_x, CUE2_y);
    disp('# of selected pairs (PRS/AOS/PC) in cell 2 = ');
    to_disp=[num2str(N_selected_PRS2),'/',num2str(N_selected_AOS2),'/',num2str(N_selected_PC2)];
    disp(to_disp);
    [rank_PRS3,N_selected_PRS3,rank_AOS3,N_selected_AOS3]=single_cell_PRS_AOS(D2D_user_list3,eNB3_x,eNB3_y,CUE3_x, CUE3_y,Max_Users3,Cell_Radius); 
    [PC_user_list3,N_selected_PC3]=single_cell_PC_Abu_channel(D2D_user_list3,eNB3_x,eNB3_y,CUE3_x, CUE3_y);
    disp('# of selected pairs (PRS/AOS/PC) in cell 3 = ');
    to_disp=[num2str(N_selected_PRS3),'/',num2str(N_selected_AOS3),'/',num2str(N_selected_PC3)];
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
%     disp('--------PRS--------');
%     disp ' '; disp('CT check, cell 1 to 2:');
    [PRS_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (PRS_user_list1, min_index12,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y); 
%     disp ' '; disp('CT check, cell 1 to 3:');                                                
    [PRS_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (PRS_user_list1_upd, min_index13,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y); 
%     disp ' '; disp('CT check, cell 2 to 1:');                                                
    [PRS_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (PRS_user_list2, min_index21,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y);
%     disp ' '; disp('CT check, cell 2 to 3:');                                                
    [PRS_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (PRS_user_list2_upd, min_index23,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
%     disp ' '; disp('CT check, cell 3 to 1:');                                                
    [PRS_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (PRS_user_list3, min_index31,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y);
%     disp ' '; disp('CT check, cell 3 to 2:');                                                
    [PRS_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (PRS_user_list3_upd, min_index32,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
    %DT forbidden region calculation
%     disp ' '; disp('DT check, cell 1 to 2:'); %FR violated -> pair will be deleted from 1st list
    [PRS_user_list1_upd,min_index12,min_index13] = forbidden_region_DUE (PRS_user_list1_upd, min_index12, ...
                                                        PRS_user_list2_upd, min_index21,DUE_SINR_min,DUE_Exp, ...
                                                                        PRS_user_list3_upd); 
%     disp ' '; disp('DT check, cell 2 to 3:'); 
    [PRS_user_list2_upd,min_index23,min_index21] = forbidden_region_DUE (PRS_user_list2_upd, min_index23, ...
                                                        PRS_user_list3_upd, min_index32,DUE_SINR_min,DUE_Exp, ...
                                                                        PRS_user_list1_upd); 
%     disp ' '; disp('DT check, cell 3 to 1:'); 
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
%     disp('--------AOS--------');
%     disp ' '; disp('CT check, cell 1 to 2:');
    [AOS_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (AOS_user_list1, min_index12,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y); 
%     disp ' '; disp('CT check, cell 1 to 3:');                                                
    [AOS_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (AOS_user_list1_upd, min_index13,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y); 
%     disp ' '; disp('CT check, cell 2 to 1:');                                                
    [AOS_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (AOS_user_list2, min_index21,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y);
%     disp ' '; disp('CT check, cell 2 to 3:');                                                
    [AOS_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (AOS_user_list2_upd, min_index23,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
%     disp ' '; disp('CT check, cell 3 to 1:');                                                
    [AOS_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (AOS_user_list3, min_index31,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y);
%     disp ' '; disp('CT check, cell 3 to 2:');                                                
    [AOS_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (AOS_user_list3_upd, min_index32,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
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
    %% -------------PC-----------------------
    %calculation of closest pair to the border
    [min_index12] = closest_pair_to_CT (PC_user_list1,CUE2_x,CUE2_y);
    [min_index13] = closest_pair_to_CT (PC_user_list1,CUE3_x,CUE3_y);
    [min_index21] = closest_pair_to_CT (PC_user_list2,CUE1_x,CUE1_y);
    [min_index23] = closest_pair_to_CT (PC_user_list2,CUE3_x,CUE3_y);
    [min_index31] = closest_pair_to_CT (PC_user_list3,CUE1_x,CUE1_y);
    [min_index32] = closest_pair_to_CT (PC_user_list3,CUE2_x,CUE2_y);

    %CT forbidden region calculation 
%     disp('--------PC--------');
%     disp ' '; disp('CT check, cell 1 to 2:');
    [PC_user_list1_upd,min_index12,min_index13] = forbidden_region_CUE (PC_user_list1, min_index12,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y); 
%     disp ' '; disp('CT check, cell 1 to 3:');                                                
    [PC_user_list1_upd,min_index13,min_index12] = forbidden_region_CUE (PC_user_list1_upd, min_index13,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y); 
%     disp ' '; disp('CT check, cell 2 to 1:');                                                
    [PC_user_list2_upd,min_index21,min_index23] = forbidden_region_CUE (PC_user_list2, min_index21,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE3_x,CUE3_y);
%     disp ' '; disp('CT check, cell 2 to 3:');                                                
    [PC_user_list2_upd,min_index23,min_index21] = forbidden_region_CUE (PC_user_list2_upd, min_index23,CUE3_x, ...
                                                        CUE3_y,eNB3_x,eNB3_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
%     disp ' '; disp('CT check, cell 3 to 1:');                                                
    [PC_user_list3_upd,min_index31,min_index32] = forbidden_region_CUE (PC_user_list3, min_index31,CUE1_x, ...
                                                        CUE1_y,eNB1_x,eNB1_y,CUE_SINR_min,CUE_Exp,CUE2_x,CUE2_y);
%     disp ' '; disp('CT check, cell 3 to 2:');                                                
    [PC_user_list3_upd,min_index32,min_index31] = forbidden_region_CUE (PC_user_list3_upd, min_index32,CUE2_x, ...
                                                        CUE2_y,eNB2_x,eNB2_y,CUE_SINR_min,CUE_Exp,CUE1_x,CUE1_y); 
    %DT forbidden region calculation
%     disp ' '; disp('DT check, cell 1 to 2:'); %FR violated -> pair will be deleted from 1st list
    [PC_user_list1_upd,min_index12,min_index13] = forbidden_region_DUE (PC_user_list1_upd, min_index12, ...
                                                        PC_user_list2_upd, min_index21,DUE_SINR_min,DUE_Exp, ...
                                                                        PC_user_list3_upd); 
%     disp ' '; disp('DT check, cell 2 to 3:'); 
    [PC_user_list2_upd,min_index23,min_index21] = forbidden_region_DUE (PC_user_list2_upd, min_index23, ...
                                                        PC_user_list3_upd, min_index32,DUE_SINR_min,DUE_Exp, ...
                                                                        PC_user_list1_upd); 
%     disp ' '; disp('DT check, cell 3 to 1:'); 
    [PC_user_list3_upd,min_index31,min_index32] = forbidden_region_DUE (PC_user_list3_upd, min_index31, ...
                                                        PC_user_list1_upd, min_index13,DUE_SINR_min,DUE_Exp, ...
                                                                        PC_user_list2_upd); 
    %% ---------------------------------------  
    random_list1=D2D_user_list1(1:15:end,:);
    random_list2=D2D_user_list2(1:15:end,:);
    random_list3=D2D_user_list3(1:15:end,:);
    random_list=[random_list1;random_list2;random_list3];    
    PRS_user_list_all_raw= [PRS_user_list1;PRS_user_list2;PRS_user_list3];
    PRS_user_list_all_new= [PRS_user_list1_upd;PRS_user_list2_upd;PRS_user_list3_upd];
    AOS_user_list_all_raw= [AOS_user_list1;AOS_user_list2;AOS_user_list3];
    AOS_user_list_all_new= [AOS_user_list1_upd;AOS_user_list2_upd;AOS_user_list3_upd];
    PC_user_list_all_raw= [PC_user_list1;PC_user_list2;PC_user_list3];
    PC_user_list_all_new= [PC_user_list1_upd;PC_user_list2_upd;PC_user_list3_upd];
    
    BS_list=[eNB1_x eNB1_y; eNB2_x eNB2_y; eNB3_x eNB3_y];
    
    CT_BS_gain1 = LTE_channel_model_urban_micro_NLOS(CUE1_x,CUE1_y,eNB1_x,eNB1_y); 
    CT_BS_gain2 = LTE_channel_model_urban_micro_NLOS(CUE2_x,CUE2_y,eNB2_x,eNB2_y);
    CT_BS_gain3 = LTE_channel_model_urban_micro_NLOS(CUE3_x,CUE3_y,eNB3_x,eNB3_y);
    
    disp('# of selected pairs after border cond. check (PRS/AOS/PC) in cell 1 = ');
    to_disp=[num2str(size(PRS_user_list1_upd,1)),'/',num2str(size(AOS_user_list1_upd,1)),'/',num2str(size(PC_user_list1_upd,1))];
    disp(to_disp);
    disp('# of selected pairs after border cond. check (PRS/AOS/PC) in cell 2 = ');
    to_disp=[num2str(size(PRS_user_list2_upd,1)),'/',num2str(size(AOS_user_list2_upd,1)),'/',num2str(size(PC_user_list2_upd,1))];
    disp(to_disp);    
    disp('# of selected pairs after border cond. check (PRS/AOS/PC) in cell 3 = ');
    to_disp=[num2str(size(PRS_user_list3_upd,1)),'/',num2str(size(AOS_user_list3_upd,1)),'/',num2str(size(PC_user_list3_upd,1))];
    disp(to_disp);
    
    %OFDMPlanning PC selection with Gain list
    [OFDMPl_PC_user_list_all,Gain_OFDMPl_PC]=OFDMPlanning_PC(D2D_user_list_all_mapped,Gain);
    [OFDMPl_PC_user_list1,OFDMPl_PC_user_list2,OFDMPl_PC_user_list3]=assign_to_cell_via_Gain(OFDMPl_PC_user_list_all,Gain_OFDMPl_PC);
    
    [Gain_OFDMPl_PC] = Gain_matrix_of_selection(D2D_user_list_all_mapped,OFDMPl_PC_user_list_all,Gain);
    [Gain_AOS] = Gain_matrix_of_selection(D2D_user_list_all_mapped,AOS_user_list_all_new,Gain);
    
    [SINR_C_OFDMPl_PC,SINR_D_OFDMPl_PC]=SINR_calc_from_Gain(Gain_OFDMPl_PC);
    [SINR_C_AOS,SINR_D_AOS]=SINR_calc_from_Gain(Gain_AOS);
    
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
    %% -----------------SINR, SE raw PC--------------------
    [SINR_C1_raw_PC,SINR_D_i1_raw_PC] = SINR_calc (1,PC_user_list1,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_raw_PC,SINR_D_i2_raw_PC] = SINR_calc (2,PC_user_list2,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_raw_PC,SINR_D_i3_raw_PC] = SINR_calc (3,PC_user_list3,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_raw,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);

    mult=1;
    for k=1:length(SINR_D_i1_raw_PC)
        mult=mult * ( 1 + SINR_D_i1_raw_PC(k) );
    end
    for k=1:length(SINR_D_i2_raw_PC)
        mult=mult * ( 1 + SINR_D_i2_raw_PC(k) );
    end
    for k=1:length(SINR_D_i3_raw_PC)
        mult=mult * ( 1 + SINR_D_i3_raw_PC(k) );
    end
    SE_raw_PC(n)=log2( (1+SINR_C1_raw_PC)*(1+SINR_C2_raw_PC)*(1+SINR_C3_raw_PC) * mult );
    %% -----------------SINR, SE multicell PC--------------------
    [SINR_C1_new_PC,SINR_D_i1_new_PC] = SINR_calc (1,PC_user_list1_upd,eNB1_x,eNB1_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C2_new_PC,SINR_D_i2_new_PC] = SINR_calc (2,PC_user_list1_upd,eNB2_x,eNB2_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);
    [SINR_C3_new_PC,SINR_D_i3_new_PC] = SINR_calc (3,PC_user_list1_upd,eNB3_x,eNB3_y,CUE1_x,CUE1_y,CUE2_x,CUE2_y,CUE3_x,CUE3_y,PC_user_list_all_new,CT_BS_gain1,CT_BS_gain2,CT_BS_gain3);

    mult=1;
    for k=1:length(SINR_D_i1_new_PC)
        mult=mult * ( 1 + SINR_D_i1_new_PC(k) );
    end
    for k=1:length(SINR_D_i2_new_PC)
        mult=mult * ( 1 + SINR_D_i2_new_PC(k) );
    end
    for k=1:length(SINR_D_i3_new_PC)
        mult=mult * ( 1 + SINR_D_i3_new_PC(k) );
    end
    SE_new_PC(n)=log2( (1+SINR_C1_new_PC)*(1+SINR_C2_new_PC)*(1+SINR_C3_new_PC) * mult );
    
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
    
    SINR_C_raw_all_PC=[SINR_C_raw_all_PC;SINR_C1_raw_PC;SINR_C2_raw_PC;SINR_C3_raw_PC];
    SINR_D_i_raw_all_PC=[SINR_D_i_raw_all_PC;SINR_D_i1_raw_PC;SINR_D_i2_raw_PC;SINR_D_i3_raw_PC];
    SINR_C_new_all_PC=[SINR_C_new_all_PC;SINR_C1_new_PC;SINR_C2_new_PC;SINR_C3_new_PC];
    SINR_D_i_new_all_PC=[SINR_D_i_new_all_PC;SINR_D_i1_new_PC;SINR_D_i2_new_PC;SINR_D_i3_new_PC];
end   

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

SINR_C_raw_all_PC_dB=pow2db(SINR_C_raw_all_PC); %PC
SINR_D_i_raw_all_PC_dB=pow2db(SINR_D_i_raw_all_PC);
SINR_C_new_all_PC_dB=pow2db(SINR_C_new_all_PC); %multiPC
SINR_D_i_new_all_PC_dB=pow2db(SINR_D_i_new_all_PC);

SE_rand=transpose(SE_rand);
SE_raw=transpose(SE_raw);
SE_new=transpose(SE_new);
SE_raw_AOS=transpose(SE_raw_AOS);
SE_new_AOS=transpose(SE_new_AOS);
SE_raw_PC=transpose(SE_raw_PC);
SE_new_PC=transpose(SE_new_PC);

%% PLOTS
% UT distribution visualization
close all;
figure
scatter(grid_Kazan(:,1),grid_Kazan(:,2),1,'MarkerEdgeAlpha',.2,'HandleVisibility','off');
% viscircles([eNB1_x eNB1_y],Cell_Radius,'Color','k','linewidth',0.5);
% viscircles([eNB2_x eNB2_y],Cell_Radius,'Color','k','linewidth',0.5);
% viscircles([eNB3_x eNB3_y],Cell_Radius,'Color','k','linewidth',0.5);
% viscircles([eNB_img_x eNB_img_y],Cell_Radius_img,'Color','k','linewidth',0.5);
axis ([-650 -50 -1025 -425]);
axis square
grid on;
hold on;
scatter([eNB1_x eNB2_x eNB3_x],[eNB1_y eNB2_y eNB3_y],'k','filled');
scatter(eNB1_x, eNB1_y,70,'r'); scatter(eNB2_x, eNB2_y,70,'b'); scatter(eNB3_x, eNB3_y,70,'m');

scatter(CT_user_list_all_mapped(:,1),CT_user_list_all_mapped(:,2),'green','filled');
scatter(CT_user_list_all_mapped(1,1), CT_user_list_all_mapped(1,2),70,'r'); 
scatter(CT_user_list_all_mapped(2,1), CT_user_list_all_mapped(2,2),70,'b'); 
scatter(CT_user_list_all_mapped(3,1), CT_user_list_all_mapped(3,2),70,'m');

scatter(D2D_user_list1(:,1),D2D_user_list1(:,2),10,'r','MarkerEdgeAlpha',.5);
scatter(D2D_user_list1(:,3),D2D_user_list1(:,4),10,'r','MarkerEdgeAlpha',.5);
scatter(D2D_user_list2(:,1),D2D_user_list2(:,2),10,'b','MarkerEdgeAlpha',.5);
scatter(D2D_user_list2(:,3),D2D_user_list2(:,4),10,'b','MarkerEdgeAlpha',.5);
scatter(D2D_user_list3(:,1),D2D_user_list3(:,2),10,'m','MarkerEdgeAlpha',.5);
scatter(D2D_user_list3(:,3),D2D_user_list3(:,4),10,'m','MarkerEdgeAlpha',.5);

scatter(OFDMPl_PC_user_list1(:,1),OFDMPl_PC_user_list1(:,2),50,'r','d');
scatter(OFDMPl_PC_user_list1(:,3),OFDMPl_PC_user_list1(:,4),50,'r','d');
scatter(OFDMPl_PC_user_list2(:,1),OFDMPl_PC_user_list2(:,2),50,'b','d');
scatter(OFDMPl_PC_user_list2(:,3),OFDMPl_PC_user_list2(:,4),50,'b','d');
scatter(OFDMPl_PC_user_list3(:,1),OFDMPl_PC_user_list3(:,2),50,'m','d');
scatter(OFDMPl_PC_user_list3(:,3),OFDMPl_PC_user_list3(:,4),50,'m','d');

scatter(AOS_user_list1_upd(:,1),AOS_user_list1_upd(:,2),30,'r','filled');
scatter(AOS_user_list1_upd(:,3),AOS_user_list1_upd(:,4),30,'r','filled');
scatter(AOS_user_list2_upd(:,1),AOS_user_list2_upd(:,2),30,'b','filled');
scatter(AOS_user_list2_upd(:,3),AOS_user_list2_upd(:,4),30,'b','filled');
scatter(AOS_user_list3_upd(:,1),AOS_user_list3_upd(:,2),30,'m','filled');
scatter(AOS_user_list3_upd(:,3),AOS_user_list3_upd(:,4),30,'m','filled');

% legend('BS','CT','init DT','DT PC','DT AOS','DT PC&AOS','DT multiPC','DT multiAOS','DT multiPC&AOS');

% %CDF DT
% figure
% cdfplot(SINR_D_i_rand_all_dB);
% xlim([0 50]);
% hold on
% % cdfplot(SINR_D_i_raw_all_dB);
% cdfplot(SINR_D_i_new_all_dB);
% % cdfplot(SINR_D_i_raw_all_AOS_dB);
% cdfplot(SINR_D_i_new_all_AOS_dB);
% % cdfplot(SINR_D_i_raw_all_PC_dB);
% cdfplot(SINR_D_i_new_all_PC_dB);
% grid on
% % legend('rand','PRS','multiPRS','AOS','multiAOS','PC','multiPC');
% legend('rand','multiPRS','multiAOS','multiPC');
% xlabel('DT SINR (dB)','FontName','Arial','FontSize',14);
% ylabel('CDF','FontName','Arial','FontSize',14);
% 
% %CDF CT
% figure
% cdfplot(SINR_C_rand_all_dB);
% xlim([0 50]);
% hold on
% % cdfplot(SINR_C_raw_all_dB);
% cdfplot(SINR_C_new_all_dB);
% % cdfplot(SINR_C_raw_all_AOS_dB);
% cdfplot(SINR_C_new_all_AOS_dB);
% % cdfplot(SINR_C_raw_all_PC_dB);
% cdfplot(SINR_C_new_all_PC_dB);
% grid on
% % legend('rand','PRS','multiPRS','AOS','multiAOS','PC','multiPC');
% legend('rand','multiPRS','multiAOS','multiPC');
% xlabel('CT SINR (dB)','FontName','Arial','FontSize',14);
% ylabel('CDF','FontName','Arial','FontSize',14);
% 
% %SE
% A=[numbers_of_pairs_rand,SE_rand]; 
% B=[numbers_of_pairs_raw,SE_raw];
% C=[numbers_of_pairs_new,SE_new];
% D=[numbers_of_pairs_raw_AOS,SE_raw_AOS];
% E=[numbers_of_pairs_new_AOS,SE_new_AOS];
% F=[numbers_of_pairs_raw_PC,SE_raw_PC];
% G=[numbers_of_pairs_new_PC,SE_new_PC];
% A=sortrows(A); %from smallest # of pairs to highest
% [Alink,aa] = findgroups(A(:,1)); %aa - different values of pairs, Alink - links from aa to A
% A1 = [ aa, splitapply(@mean,A(:,2),Alink)]; %find mean value for each # of pairs
% B=sortrows(B);
% [Blink,bb] = findgroups(B(:,1));
% B1 = [ bb, splitapply(@mean,B(:,2),Blink)];
% C= sortrows(C);
% [Clink,cc] = findgroups(C(:,1));
% C1 = [ cc, splitapply(@mean,C(:,2),Clink)];
% D= sortrows(D);
% [Dlink,dd] = findgroups(D(:,1));
% D1 = [ dd, splitapply(@mean,D(:,2),Dlink)];
% E= sortrows(E);
% [Elink,ee] = findgroups(E(:,1));
% E1 = [ ee, splitapply(@mean,E(:,2),Elink)];
% F= sortrows(F);
% [Flink,ff] = findgroups(F(:,1));
% F1 = [ ff, splitapply(@mean,F(:,2),Flink)];
% G= sortrows(G);
% [Glink,gg] = findgroups(G(:,1));
% G1 = [ gg, splitapply(@mean,G(:,2),Glink)];
% 
% figure
% plot(A1(:,1),A1(:,2),'b-o','linewidth',2.5); %b r y m g c 
% hold on
% % plot(B1(:,1),B1(:,2),'r-o','linewidth',2.5);
% plot(C1(:,1),C1(:,2),'r-o','linewidth',2.5);
% % plot(D1(:,1),D1(:,2),'m-o','linewidth',2.5);
% plot(E1(:,1),E1(:,2),'y-o','linewidth',2.5);
% % plot(F1(:,1),F1(:,2),'c-o','linewidth',2.5);
% plot(G1(:,1),G1(:,2),'m-o','linewidth',2.5);
% grid on
% % legend('rand','PRS','multiPRS','AOS','multiAOS','PC','multiPC');
% legend('rand','multiPRS','multiAOS','multiPC');
% xlabel('Number of D2D pairs','FontName','Arial','FontSize',14);
% ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);

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
%             disp('list is empty');
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
                                        PRS_user_list2, min_index2,SINR,Exp,PRS_user_list_other) 
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
        r1=dist(PRS_user_list1(min_index1,1),PRS_user_list1(min_index1,2), ... %pair1 Tx - pair2 Rx ;was pair1 Rx - pair2 Rx in ARP
            PRS_user_list2(min_index2,3),PRS_user_list2(min_index2,4)); 
        R= (SINR ^(2/Exp)) * l * r0 / r1;
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
    SINR_th_BS_W=db2pow(SINR_th_BS); %10dB = 10 Watts
    SINR_th_DT_W=db2pow(SINR_th_DT); %15dB = 31.6228 Watts
    
    Pd = 1*10^(-3); %W
    Pc = 50*10^(-3); %W
    Bandwidth_kHz = 180; 
    Noise_dBm =	-174; 
    Noise_Total_Watts =	Bandwidth_kHz*10^3*(db2pow(Noise_dBm-30)); 
    while 1
        %CT
        DT_rx_SINR=[];
        sum1=0; %sum over DTs interference towards BS      
        for i=1:size(D2D_user_list,1)
            sum1= sum1 + (Pd * LTE_channel_model_urban_micro_NLOS(D2D_user_list(i,1),D2D_user_list(i,2),eNB_x,eNB_y));
            % LTE_channel_model_urban_micro_NLOS is gain, gain=1/loss, so there is *    
        end
        BS_rx_SINR = Pc * LTE_channel_model_urban_micro_NLOS(CUE_x,CUE_y,eNB_x,eNB_y) / (sum1 + Noise_Total_Watts);
        
        %DT
        for i=1:size(D2D_user_list,1)
            sum2=0; %sum over DTs interference towards i-th DTrx    
            for j=1:size(D2D_user_list,1)
                if (j~=i)
                    sum2= sum2 + Pd * LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(j,1),D2D_user_list(j,2),D2D_user_list(i,3),D2D_user_list(i,4));           
                end
            end
            DT_rx_SINR(i) = Pd * LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(i,1),D2D_user_list(i,2),D2D_user_list(i,3),D2D_user_list(i,4)) / ...
                        ( Pc * LTE_channel_model_indoor_hotspot_NLOS(CUE_x,CUE_y,D2D_user_list(i,3),D2D_user_list(i,4)) + sum2 + Noise_Total_Watts);  %was urban
        end 
        
        DT_rx_SINR=transpose(DT_rx_SINR);
        minimum = min(DT_rx_SINR); 
        row=find(DT_rx_SINR==minimum);
        if ( (isempty(D2D_user_list))  || ( (BS_rx_SINR>SINR_th_BS_W) &&  (DT_rx_SINR(row)>SINR_th_DT_W) ) ) %
            N_selected_PC=size(D2D_user_list,1);
            break
        else
            D2D_user_list(row,:)=[];
        end
    end
end

function [D2D_user_list_all,Gain] = OFDMPlanning_PC(D2D_user_list_all,Gain)
    %Gain should be in times
    SINR_th_CT_dB=-30; %dB
    SINR_th_DT_dB=-40; %dB
    % SINR_th_CT=db2pow(SINR_th_CT_dB); %Watts; 6dB = 3.9811 
    % SINR_th_DT=db2pow(SINR_th_DT_dB); %Watts; 15dB = 31.6228 
    while 1
        [SINR_C_dB,SINR_D_dB] = SINR_calc_from_Gain (Gain);
        SINR_D_dB_Tx=SINR_D_dB(1:2:end,:);
        min_row_DT=find(SINR_D_dB_Tx==min(SINR_D_dB_Tx));
        min_row_CT=find(SINR_C_dB==min(SINR_C_dB));
        if ( (isempty(D2D_user_list_all))  || ( (SINR_D_dB_Tx(min_row_DT)>SINR_th_DT_dB) &&  (SINR_C_dB(min_row_CT)>SINR_th_CT_dB) ) )   
            break
        else
            D2D_user_list_all(min_row_DT,:)=[];
            SINR_D_dB_Tx(min_row_DT,:)=[];
            Gain(2*min_row_DT-1+3,:)=[]; %remove Tx
            Gain(:,2*min_row_DT-1+3+3)=[];
            Gain(2*min_row_DT-1+3,:)=[]; %remove Rx, here Rx comes to the Tx's place
            Gain(:,2*min_row_DT-1+3+3)=[]; 
        end
    end
end

function [SINR_C_dB,SINR_D_dB] = SINR_calc_from_Gain (Gain)
    number_of_DTs=size(Gain,1)-6;
    if number_of_DTs==0
        SINR_D=[];SINR_D_dB=[];
    end
    CT_BS_gain_all = [Gain(1,1) Gain(2,2) Gain(3,3)];
    DT_BS_gain_all=Gain(4:number_of_DTs+3,1:3);
    CT_DT_gain_all=Gain(4:number_of_DTs+3,4:6);
    DT_DT_gain_all=Gain(4:number_of_DTs+3,7:number_of_DTs+6); %incl. 0 and DT_Pair_gain
    DT_Pair_gain=[];
    for i=1:number_of_DTs
        switch rem(i,2) %remainder
            case 1 %odd
                DT_Pair_gain(i)=Gain(i+3,i+7);
            case 0 %even
                DT_Pair_gain(i)=Gain(i+3,i+5);
       end
    end
    DT_Pair_gain=transpose(DT_Pair_gain);    
    Bandwidth_kHz = 180; 
    Noise_dBm =	-174; 
    Noise_Total_Watts =	Bandwidth_kHz*10^3*(db2pow(Noise_dBm-30)); 
    DT_Tx_Power = 1*10^(-3); %DT_Tx_Power=1*10^(-3) W; 10*log10(1) dBm
    CT_Tx_Power = 50*10^(-3); 
    
    %CT
    for numb_cell=1:3
        switch numb_cell
            case 1
                CT_BS_gain=CT_BS_gain_all(1);
                DT_BS_gain=DT_BS_gain_all(:,1);
                CT_BS_gain_other_1 = Gain(2,1); %eNB1 to CT2 
                CT_BS_gain_other_2 = Gain(3,1); %eNB1 to CT3 
            case 2
                CT_BS_gain=CT_BS_gain_all(2);
                DT_BS_gain=DT_BS_gain_all(:,2);
                CT_BS_gain_other_1 = Gain(1,2);  
                CT_BS_gain_other_2 = Gain(3,2); 
            case 3 
                CT_BS_gain=CT_BS_gain_all(3); 
                DT_BS_gain=DT_BS_gain_all(:,3);
                CT_BS_gain_other_1 = Gain(1,3);
                CT_BS_gain_other_2 = Gain(2,3);
        end                
        int_from_CTs_to_CT=0; %interference
        int_from_CTs_to_CT=CT_Tx_Power*CT_BS_gain_other_1+CT_Tx_Power*CT_BS_gain_other_2; 
        int_from_DTs_to_CT=0;
        for i=1:size(DT_BS_gain,1)
            int_from_DTs_to_CT=int_from_DTs_to_CT+ (DT_Tx_Power*DT_BS_gain(i));
        end    
        SINR_C(numb_cell) = (CT_Tx_Power*CT_BS_gain) / (Noise_Total_Watts + int_from_DTs_to_CT +int_from_CTs_to_CT);%times        
        SINR_C_dB(numb_cell)=pow2db(SINR_C(numb_cell));
    end
    SINR_C=transpose(SINR_C);%times
    SINR_C_dB=transpose(SINR_C_dB);
    
    %DT
    for i=1:size(DT_DT_gain_all,1)
        int_from_CTs_to_DT=0;
        int_from_CTs_to_DT=CT_Tx_Power*CT_DT_gain_all(i,1)+CT_Tx_Power*CT_DT_gain_all(i,2)+CT_Tx_Power*CT_DT_gain_all(i,3); 
        int_from_DTs_to_DT=0; 
        for j=1:(size(DT_DT_gain_all,2))
            if i~=j
                int_from_DTs_to_DT=int_from_DTs_to_DT+(DT_Tx_Power*DT_DT_gain_all(i,j)); 
            end    
        end
        SINR_D(i) = (DT_Tx_Power*DT_Pair_gain(i))/( int_from_DTs_to_DT +int_from_CTs_to_DT + Noise_Total_Watts);
        SINR_D_dB(i) = pow2db(SINR_D(i));
    end
    SINR_D=transpose(SINR_D);%times
    SINR_D_dB=transpose(SINR_D_dB);
end

function [D2D_user_list1,D2D_user_list2,D2D_user_list3]= assign_to_cell_via_Gain(D2D_user_list_all_mapped,Gain)
    Gain_DTs=Gain(4:end-3,1:3); %only DTs needed
    D2D_user_list1=[];D2D_user_list2=[];D2D_user_list3=[];
    k=1;
    for i=1:size(D2D_user_list_all_mapped,1)
        [val,ind]=max(Gain_DTs(k,:));
        switch ind %assign pair to BS according to DT Tx
            case 1
                D2D_user_list1=[D2D_user_list1 ; D2D_user_list_all_mapped(i,:)];
                k=k+2; %to take only odd rows of Gain_DTs (DT Tx)
            case 2
                D2D_user_list2=[D2D_user_list2 ; D2D_user_list_all_mapped(i,:)];
                k=k+2;
            case 3
                D2D_user_list3=[D2D_user_list3 ; D2D_user_list_all_mapped(i,:)];
                k=k+2;
        end
    end    
end

function [Gain] = Gain_matrix_of_selection(D2D_user_list_all_mapped,selection_user_list,Gain)
    [~,numbers_of_selected_pairs,~] = intersect(D2D_user_list_all_mapped,selection_user_list,'rows');
    numbers_of_selected_pairs=sort(numbers_of_selected_pairs,'desc');  
    k=1;
    for i=size(D2D_user_list_all_mapped,1):-1:1
        if i~=numbers_of_selected_pairs(k)
            Gain(3+2*i-1,:)=[]; %remove Tx
            Gain(:,3+2*i-1+3)=[];
            Gain(3+2*i-1,:)=[]; %remove Rx, here Rx comes to the Tx's place
            Gain(:,3+2*i-1+3)=[]; 
        else 
            if k~=length(numbers_of_selected_pairs)
                k=k+1;
            end
        end   
    end
end