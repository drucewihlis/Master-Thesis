function [rank_PRS,N_selected_PRS,rank,N_selected_AOS]=single_cell_PRS_AOS(D2D_user_list, eNB_x,eNB_y,CUE_x_tx,CUE_y_tx,Max_Users,Cell_Radius)

% clc;clear all;close all;
%% ParamIters Definitions as per table 1
Noise_dBm =	-174; % AWG Noise for eNB, CUE & DUE in dBm/Hz
CT_SNR_Target_dB = 12; %dB
DT_SNR_Target_dB = 25; %dB
CUE_SINR_min_dB	= 6; %dB
DUE_SINR_min_dB	= 15; %dB
UE_P_Tx_Max_dBm = 23; % Maximum UE transmit power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 72
UE_P_Tx_Min_dBm = -40; % Minimum UE transmit power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 115
CT_P_Rx_min_dBm = -98; % Minimum CT receiver power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 222
DT_P_Rx_min_dBm = -78; % Minimum D2D ProSe receiver power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 222
Carrier_Frequency =	2.6; % GHz
Bandwidth_kHz = 180; % Bandwidth per resource block (RB) in LTE
CT_Tx_Power_UIP_dBm = 23; % set CT power at which the value of k saturates the SE

%% Conversion of above paramIters to ratio or Watts values
Noise_Total_Watts =	Bandwidth_kHz*10^3*(db2pow(Noise_dBm-30)); % Total AWG Noise for eNB, CUE & DUE in Watts
Noise_Total_dBm = pow2db(Noise_Total_Watts)+30; % Total AWG Noise for eNB, CUE & DUE in dBm
CT_SNR_Target = db2pow(CT_SNR_Target_dB); %Real ratio value
DT_SNR_Target = db2pow(DT_SNR_Target_dB); %Real ratio value
CUE_SINR_min = db2pow(CUE_SINR_min_dB); % Real ratio value
DUE_SINR_min = db2pow(DUE_SINR_min_dB); % Real ratio value
UE_P_Tx_Max = db2pow(UE_P_Tx_Max_dBm-30); % Maximum UE transmit power in Watts
UE_P_Tx_Min = db2pow(UE_P_Tx_Min_dBm-30); % Minimum UE transmit power in Watts
CT_P_Rx_min = db2pow(CT_P_Rx_min_dBm-30); % Minimum CT receiver power in Watts
DT_P_Rx_min = db2pow(DT_P_Rx_min_dBm-30); % Minimum D2D ProSe receiver power in Watts

CUE_Exp = 3.67; % Cellular Pathloss Exponent - Alpha_zero
DUE_Exp = 4.33; % Devices Pathloss Exponent -Alpha_d
CUE_Prop_Const = (10^-2.27)*(Carrier_Frequency^-2.6); % Cellular Propagation Constant - C_zero
DUE_Prop_Const = (10^-1.15)*(Carrier_Frequency^-2); % Devices Propagation Constant - C_d
CUE_Prop_Const_dB = pow2db(CUE_Prop_Const); % Cellular Propagation Constant - C_zero in dB
DUE_Prop_Const_dB = pow2db(DUE_Prop_Const); % Devices Propagation Constant - C_d in dB

% Cell_Radius = 200; % Cell Radius in meters
UE_Dist_Min = 10; % Minimum distance of any UE (i.e. CUE or DUE) from either the BS or another UE (i.e. CUE or DUE)
CUE_DUE_rx_Dist_Max = 1 * Cell_Radius; % Maximum distance between CUE & receivers of other DUEs
DUE_DUE_rx_Dist_Max = 1 * Cell_Radius; % Maximum distance between DUE & receivers of other DUEs

D2D_Sep_Max = 0.1*Cell_Radius;
% Max_Users = 100;
No_Runs = 1; %initially 10

CT_Int_Power_UIP = zeros(Max_Users,No_Runs);
DT_SINR_3_UIP_dB = zeros (Max_Users,No_Runs);

for Iter = 1:No_Runs
    %% Placement of Users
   % eNB_x = Cell_Radius; 
   % eNB_y = Cell_Radius;
%     locUE(Iter) = UE_Dist_Min + (Cell_Radius - UE_Dist_Min)*sqrt(rand(1,1));
    % Generate the random angle Theta of the points:
%     theta_= 2*pi*rand(1,1);
%     CUE_x_tx = locUE(Iter)*cos(theta_) + eNB_x ;
%     CUE_y_tx = locUE(Iter)*sin(theta_) + eNB_y ;
    R_CT_eNB = pdist([CUE_x_tx CUE_y_tx;eNB_x eNB_y]);
    
    %% Create forbidden region
%     Rforb_c= CT_P_Rx_min^(2/CUE_Exp)*l*r0/r1;
 
    %% Threshold Distances to Control the selection Criteria
    R_CT_DT_Reuse = (CT_SNR_Target*DUE_SINR_min/(DT_SNR_Target-DUE_SINR_min))^(1/DUE_Exp)*R_CT_eNB;
    d_DT_DT_Reuse = (DT_SNR_Target*DUE_SINR_min/(DT_SNR_Target-DUE_SINR_min))^(1/DUE_Exp)*D2D_Sep_Max;
    R_DT_eNB_Reuse = (DT_SNR_Target*CUE_SINR_min/(CT_SNR_Target-CUE_SINR_min))^(1/DUE_Exp)*D2D_Sep_Max;
 
    %% Separation Distances between all communicating entities (CT, eNB and DTs) 
%     D2D_user_list =
%     LTE_UE_uniform_distribution(eNB_x,eNB_y,Cell_Radius,D2D_Sep_Max, Max_Users); 
    %is gererated in multi_cell_PRS
    Link_list = D2D_user_list;
    Link_list(Max_Users+1,:) = [CUE_x_tx, CUE_y_tx, eNB_x, eNB_y];
    for j=1:Max_Users+1
        for k=1:Max_Users+1
            Distance(j,k) = pdist([Link_list(j,1) Link_list(j,2);Link_list(k,3) Link_list(k,4)]);
        end
    end
    %% AOS-Power allocation according to SNR Target
    CT_BS_gain = LTE_channel_model_urban_micro_NLOS(CUE_x_tx,CUE_y_tx,eNB_x,eNB_y);
    CT_Tx_Power_AOS = CT_SNR_Target*Noise_Total_Watts/CT_BS_gain;
    CT_Rx_Power_AOS = CT_Tx_Power_AOS * CT_BS_gain;
 
    %% Link Gains & Interference Power Computations
    for jj=1:Max_Users
        DT_Pair_gain(jj) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(jj,1),D2D_user_list(jj,2),D2D_user_list(jj,3),D2D_user_list(jj,4));
        P_DT_AOS(jj) = (DT_SNR_Target*Noise_Total_Watts)/DT_Pair_gain(jj); % Transmit power set based on target SNR at the DT receiver
    end
    for jj = 1:Max_Users
        % DT received power
        DT_rx_Power_AOS(jj) = P_DT_AOS(jj) * DT_Pair_gain(jj);
        % CT interference power from the DT
        DT_eNB_gain_AOS(jj) = LTE_channel_model_urban_micro_NLOS(D2D_user_list(jj,1),D2D_user_list(jj,2),eNB_x,eNB_y);
        %DT_eNB_gain_AOS(jj) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(jj,1),D2D_user_list(jj,2),eNB_x,eNB_y);
        CT_Int_Power_from_DT_AOS(jj) = P_DT_AOS(jj) * DT_eNB_gain_AOS(jj);
        % DT intereference power from CT
        DT_CT_gain(jj) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(jj,3),D2D_user_list(jj,4),CUE_x_tx,CUE_y_tx);
        DT_Int_Power_from_CT_AOS(jj) = CT_Tx_Power_AOS * DT_CT_gain(jj);
    end
    % DT intereference power from other DTs excluding CT
    for jj = 1:Max_Users
        for tt = 1:Max_Users
            if tt ~= jj
                DT_DT_gain(jj,tt) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(jj,3),D2D_user_list(jj,4),D2D_user_list(tt,1),D2D_user_list(tt,2));
                DT_Int_Power_from_DT_AOS(jj,tt) = P_DT_AOS(tt)*DT_DT_gain(jj,tt);
            end
        end
    end
%%-------------- End of common part ----------------------------- 
    %% ---------- AOS ----------------
    % Rank DTs according to priority to select the first DT sharing CT RB
    SE_Initial = log2(1+CT_SNR_Target);
    N_selected_AOS(Iter) = 0;
    CT_Final_SINR_1_AOS(Iter) = 0;
    rank = zeros(Max_Users,1);
    DT_Ranklist_AOS = zeros (Max_Users,1);
    ds=zeros(Max_Users,Max_Users);
    w1 = 1; w2 = 1; w3 = 1; d1 = zeros(Max_Users);
    d1 = w1*Distance(Max_Users+1,1:Max_Users) + w2*Distance(1:Max_Users,Max_Users+1)';    
    v=0;am=0;
    for ids=1:Max_Users
        b1 = Distance(Max_Users+1,ids) ; 
        b2 = Distance(ids,Max_Users+1) ;
        if and (b1 > R_CT_DT_Reuse , b2 > R_DT_eNB_Reuse)
            if d1(ids)> v
                v= d1(ids);
                am =ids; 
            end
        else
            rank(ids)=-1;
        end
    end
    ranked = 1;
    rank(am) = ranked;
    DT_Ranklist_AOS(ranked) = am;
    N_selected_AOS(Iter) = N_selected_AOS(Iter) + 1;
    %% SINR for CT & first selected DT AOS Approach
    CT_Total_Int_Power_AOS(N_selected_AOS(Iter),Iter) = CT_Int_Power_from_DT_AOS(am);
    CT_SINR_1_AOS(N_selected_AOS(Iter),Iter) = CT_Rx_Power_AOS/(CT_Total_Int_Power_AOS(N_selected_AOS(Iter),Iter) + Noise_Total_Watts);
    DT_SINR_1_AOS(am) = DT_rx_Power_AOS(am)/(DT_Int_Power_from_CT_AOS(am)+Noise_Total_Watts);
    SE_1_AOS(N_selected_AOS(Iter),Iter)= log2(1+CT_SINR_1_AOS(N_selected_AOS(Iter),Iter)) + log2(1+DT_SINR_1_AOS(am));
    SE_prev = SE_Initial;    
    ds(:,ranked)=d1;
    %% Selection of subsquent DT according to priority
    while and (SE_1_AOS(N_selected_AOS(Iter),Iter) > SE_prev, ranked < Max_Users)
        ranked= ranked+1;
        d2 = zeros(Max_Users,1);
        candidate=0;
        for ids=1:Max_Users
            if rank(ids)==0 %check unranked DTs
                past = 1;
                d2(ids)= d1(ids);
                C1 = Distance(Max_Users+1,ids) < (N_selected_AOS(Iter)+1)^(1/DUE_Exp)*R_CT_DT_Reuse; 
                C2 = Distance(ids,Max_Users+1) < (N_selected_AOS(Iter)+1)^(1/DUE_Exp)*R_DT_eNB_Reuse; 
                for idrs=1:Max_Users
                    if rank(idrs) > 0
                        C3 = Distance(ids,idrs) < (N_selected_AOS(Iter)+1)^(1/DUE_Exp)*d_DT_DT_Reuse;
                        if or( or (C1, C2), C3)
                            d2(ids)=0;
                            rank(ids)=-1;
                            past = 0;
                            break
                        end
                        d2(ids)=d2(ids) + w3 * Distance(ids,idrs)';
                    end
                end
                if past == 1
                    candidate = 1;
                end  
            end
        end
        %get the best next Nth DT to share CT RB     
        if(candidate~=0)
            v2=0;am2=0;
            for ids=1:Max_Users
                if d2(ids)> v2
                    v2=d2(ids);
                    am2 =ids; 
                end
            end
            if am2 ~=0
                rank(am2) = ranked;
                DT_Ranklist_AOS(ranked) = am2;
                N_selected_AOS(Iter) = N_selected_AOS(Iter) + 1;
                %% SINR for CT & subsquent selected DTs AOS Approach
                CT_Total_Int_Power_AOS(N_selected_AOS(Iter),Iter) = CT_Total_Int_Power_AOS(N_selected_AOS(Iter)-1,Iter) + CT_Int_Power_from_DT_AOS(am2);
                CT_SINR_1_AOS(N_selected_AOS(Iter),Iter) = CT_Rx_Power_AOS/(CT_Total_Int_Power_AOS(N_selected_AOS(Iter),Iter) + Noise_Total_Watts);
                %% Compute SINRs everytime new DT is opportunistically selected
                DT_Total_Int_Power_AOS = zeros(Max_Users,1);
                for ii=1:N_selected_AOS(Iter)
                    rr = DT_Ranklist_AOS(ii);
                    DT_Total_Int_Power_AOS(rr)= DT_Int_Power_from_CT_AOS(rr);
                    for jj = 1:N_selected_AOS(Iter)
                        tt = DT_Ranklist_AOS(jj);
                        if (jj~=ii)    
                            DT_Total_Int_Power_AOS(rr)= DT_Total_Int_Power_AOS(rr)+ DT_Int_Power_from_DT_AOS(rr,tt);
                        end
                    end
                end
                %% Computation of SE everytime a new DT is opportunistically selected
                SE_prev = SE_1_AOS(N_selected_AOS(Iter)-1,Iter);
                SE_1_AOS(N_selected_AOS(Iter),Iter)= log2(1+CT_SINR_1_AOS(N_selected_AOS(Iter),Iter));
                DT_SINR_1_AOS = zeros(Max_Users,1);
                for ii=1:Max_Users
                    if rank(ii) > 0    
                        DT_SINR_1_AOS(ii) = DT_rx_Power_AOS(ii)/(DT_Total_Int_Power_AOS(ii)+Noise_Total_Watts);
                        DT_SINR_1_AOS_dB(ii,Iter) = pow2db(DT_SINR_1_AOS(ii));
                        SE_1_AOS(N_selected_AOS(Iter),Iter)=SE_1_AOS(N_selected_AOS(Iter),Iter)+log2(1+DT_SINR_1_AOS(ii));
                    end
                end
                ds(:,ranked)=d2;
            end        
        end
    end
    %    decrement N_selected_AOS
    SE_final_AOS = SE_prev;
    CT_Final_SINR_1_AOS(Iter) = CT_SINR_1_AOS(N_selected_AOS(Iter),Iter);
    CT_Final_SINR_1_AOS_dB(Iter) = pow2db(CT_Final_SINR_1_AOS(Iter));

   
    %% -------- Random DT Selection PRS Approach ----------------
    % Initialization of some value
    N_selected_PRS(Iter)=0;
    rank_PRS = zeros(Max_Users,1);
    CT_Final_SINR_2_PRS(Iter) = 0;
    DT_Ranklist_PRS = zeros (Max_Users,1); 
    ranked_PRS = 0;
    SE_prev_PRS = SE_Initial;    
    %% Random selection of DT
    % Select randomly any of the DTs that don't satisfy C1, C2, & C3 (with N_selected_PRS==0)
    Random_DT_PRS = randperm(Max_Users);
    for ranked_PRS = 1:Max_Users
        candidate_found_PRS = 0;
        am_PRS=0;
        for ids_PRS=1:Max_Users
            past_PRS =0;
            cc=Random_DT_PRS(ids_PRS);
            if rank_PRS(cc)==0 %check unranked DTs
                C1 = Distance(Max_Users+1,cc) < (N_selected_PRS(Iter)+1)^(1/DUE_Exp)*R_CT_DT_Reuse; 
                C2 = Distance(cc,Max_Users+1) < (N_selected_PRS(Iter)+1)^(1/DUE_Exp)*R_DT_eNB_Reuse; 
                past_PRS = 1;
                for idrs_PRS=1:Max_Users
                    if and(rank_PRS(idrs_PRS) > 0, idrs_PRS ~= cc)
                        C3 = Distance(cc,idrs_PRS) < (N_selected_PRS(Iter)+1)^(1/DUE_Exp)*d_DT_DT_Reuse;
                        if or( or (C1, C2), C3)
                            rank_PRS(cc)=-1;
                            past_PRS = 0;
                            break
                        end
                    end
                end
            end   
            if past_PRS == 1
                candidate_found_PRS = 1;
                am_PRS = cc;
                break
            end
        end
        if(candidate_found_PRS == 1)
            rank_PRS(am_PRS) = ranked_PRS;
            DT_Ranklist_PRS(ranked_PRS) = am_PRS;
            N_selected_PRS(Iter) = N_selected_PRS(Iter) + 1;
            %% Compute SINRs every time an Nth DT is randomly selected
            % SINR for CT & subsquent selected DTs PRS Approach
            if N_selected_PRS(Iter) ==1
                CT_Total_Int_Power_PRS(N_selected_PRS(Iter),Iter) = CT_Int_Power_from_DT_AOS(am_PRS);
            else
                CT_Total_Int_Power_PRS(N_selected_PRS(Iter),Iter) = CT_Total_Int_Power_PRS(N_selected_PRS(Iter)-1,Iter) + CT_Int_Power_from_DT_AOS(am_PRS);
            end
            CT_SINR_2_PRS(N_selected_PRS(Iter),Iter) = CT_Rx_Power_AOS/(CT_Total_Int_Power_PRS(N_selected_PRS(Iter),Iter) + Noise_Total_Watts);
            %% Compute DTs' SINRs everytime new DT is randomly selected
            DT_Total_Int_Power_PRS = zeros(Max_Users,1);
            for ii=1:N_selected_PRS(Iter)
                rr = DT_Ranklist_PRS(ii);
                DT_Total_Int_Power_PRS(rr)= DT_Int_Power_from_CT_AOS(rr); % this power is picked from same interference matrix
                if N_selected_PRS(Iter) > 1
                    for jj = 1:N_selected_PRS(Iter)
                        tt = DT_Ranklist_PRS(jj);
                        if (jj~=ii)    
                            DT_Total_Int_Power_PRS(rr)= DT_Total_Int_Power_PRS(rr)+ DT_Int_Power_from_DT_AOS(rr,tt);
                        end
                    end
                end
            end
            %% Computation of SE everytime a new DT is randomly selected
            if N_selected_PRS(Iter) > 1
                SE_prev_PRS = SE_2_PRS(N_selected_PRS(Iter)-1,Iter);
            end
            SE_2_PRS(N_selected_PRS(Iter),Iter)= log2(1+CT_SINR_2_PRS(N_selected_PRS(Iter),Iter));
            DT_SINR_2_PRS = zeros(Max_Users,1);
            for ii=1:Max_Users
                if rank_PRS(ii) > 0
                    DT_SINR_2_PRS(ii) = DT_rx_Power_AOS(ii)/(DT_Total_Int_Power_PRS(ii)+Noise_Total_Watts);
                    DT_SINR_2_PRS_dB(ii,Iter) = pow2db(DT_SINR_2_PRS(ii));
                    SE_2_PRS(N_selected_PRS(Iter),Iter)=SE_2_PRS(N_selected_PRS(Iter),Iter)+log2(1+DT_SINR_2_PRS(ii));
                end
            end
        end
        if (SE_2_PRS(N_selected_PRS(Iter),Iter) < SE_prev_PRS)
            break
        end
    end         
    %    decrement N_Selected_PRS & return the final N_selected_PRS
    SE_final_PRS = SE_prev_PRS;
    CT_Final_SINR_2_PRS(Iter) = CT_SINR_2_PRS(N_selected_PRS(Iter),Iter);
    CT_Final_SINR_2_PRS_dB(Iter) = pow2db(CT_Final_SINR_2_PRS(Iter));


    %% -------- UIP scheme approach using above N_selected_AOS ----------
    %deleted%
end
%% Computation of mean SE for each row as DTs are added
%AOS, UIP deleted
size_2_PRS = size(SE_2_PRS);
CleanSE_2_PRS = zeros(size_2_PRS);
[Val_2_PRS,Loc_2_PRS] = max(SE_2_PRS);
for tc_PRS = 1:No_Runs
    CleanSE_2_PRS(Loc_2_PRS(tc_PRS),tc_PRS)=Val_2_PRS(tc_PRS);
end
SE_cleanmean_2_PRS = sum(CleanSE_2_PRS,2) ./ sum(CleanSE_2_PRS~=0,2); % comma 2 indicates average over rows otherwise use comma 1 to take average over columns
SE_mean_2_PRS = sum(SE_2_PRS,2) ./ sum(SE_2_PRS~=0,2); % comma 2 indicates average over rows otherwise use comma 1 to take average over columns
Size_meanSE_2_PRS = size(SE_mean_2_PRS,1);
N_Selected_max = min([Size_meanSE_2_PRS]);
N_DTs = 1:N_Selected_max;
%% Plot of SE vs N_DTs for AOS, PRS & UIP
% figure
% hold on
% plot(N_DTs,SE_cleanmean_2_PRS(N_DTs),'b-o','linewidth',2.5);
% set(gca,'fontsize',14)
% grid on 
% xlabel('Number of D2D Pairs','FontName','Arial','FontSize',14);
% ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);
% xlim([0 15]);
% %title('SE vs No. D2D Pairs');
% 
% %% Plot of second SE vs N_DTs for AOS, PRS & UIP
% figure
% hold on
% plot(N_DTs,SE_mean_2_PRS(N_DTs),'b-o','linewidth',2.5);
% set(gca,'fontsize',14)
% grid on 
% xlabel('Number of D2D Pairs','FontName','Arial','FontSize',14);
% ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);
% xlim([0 15]);
% title('Second SE vs No. D2D Pairs');
% 
% %% Plot of CDF for realised Spectral Efficiency over multiple runs
% SE_2_PRS(SE_2_PRS==0) = []; % Removal of all zero values 
% figure
% hold on
% grid on
% A8=cdfplot(SE_2_PRS(:));
% A8.LineWidth = 2.5; 
% set(A8,'Color','b','LineStyle','-');
% set(gca,'fontsize',14);
% xlabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);
% ylabel('CDF','FontName','Arial','FontSize',14);
% %xlim([0 50]);
% hold off
% 
% %% Plot of CDF for realised CT SINRs over multiple runs
% figure
% hold on
% grid on
% A2=cdfplot(CT_Final_SINR_2_PRS_dB(:));
% A2.LineWidth = 2.5; 
% set(A2,'Color','b','LineStyle','-');
% set(gca,'fontsize',14);
% xlabel('CT SINR (dB)','FontName','Arial','FontSize',14);
% ylabel('CDF','FontName','Arial','FontSize',14);
% xlim([0 15]);
% hold off
% 
% %% Plot of CDF for realised DTs SINRs over multiple runs
% % DT_SINR_1_AOS_dB(DT_SINR_1_AOS_dB==0) = []; % Removal of all zero values 
% DT_SINR_2_PRS_dB(DT_SINR_2_PRS_dB==0) = []; % Removal of all zero values 
% % DT_SINR_3_UIP_dB(DT_SINR_3_UIP_dB==0) = []; % Removal of all zero values 
% figure
% hold on
% grid on
% A5=cdfplot(DT_SINR_2_PRS_dB(:));
% A5.LineWidth = 2.5; 
% set(A5,'Color','b','LineStyle','-');
% line([15 15],ylim,'Color','k','LineStyle','-.');
% set(gca,'fontsize',14);
% xlabel('DTs SINR (dB)','FontName','Arial','FontSize',14);
% ylabel('CDF','FontName','Arial','FontSize',14);
% xlim([0 50]);
% hold off
% 
% %% Plot of CDF for selected number of D2D pairs over multiple runs
% figure
% A10.LineWidth = 2.5; 
% hold on
% grid on
% A11=cdfplot(N_selected_PRS(:));
% A11.LineWidth = 2.5; 
% set(A11,'Color','b','LineStyle','-');
% set(gca,'fontsize',14);
% xlabel('Number of D2D Pairs','FontName','Arial','FontSize',14);
% ylabel('CDF','FontName','Arial','FontSize',14);
% xlim([0 15]);
% hold off

end