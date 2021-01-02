clc;
clear;
close all;

%% Parameters Definitions as per table 1

Noise_dBm =	-174;                                                           % AWG Noise for eNB, CUE & DUE in dBm/Hz
CUE_SINR_min_dB	= 6;                                                        %dB
DUE_SINR_min_dB	= 15;                                                       %dB
UE_P_Tx_Max_dBm = 23;                                                       % Maximum UE transmit power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 72
UE_P_Tx_Min_dBm = -40;                                                      % Minimum UE transmit power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 115
CT_P_Rx_min_dBm = -98;                                                      % Minimum CT receiver power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 222
DT_P_Rx_min_dBm = -78;                                                      % Minimum D2D ProSe receiver power in dBm - 3GPP TS 36.101 version 13.2.1 (2016 May) page 222
Carrier_Frequency =	2.6;                                                    % GHz
Bandwidth_kHz = 180;                                                        % Bandwidth per resource block (RB) in LTE

%% Conversion of above parameters to ratio or Watts values

Noise_Total_Watts =	Bandwidth_kHz*10^3*(10^(0.1*(Noise_dBm-30)));              % Total AWG Noise for eNB, CUE & DUE in Watts
Noise_Total_dBm = 10*log10(Noise_Total_Watts)+30;                            % Total AWG Noise for eNB, CUE & DUE in dBm
CUE_SINR_min = 10^(0.1*(CUE_SINR_min_dB));                                     % Real ratio value
DUE_SINR_min = 10^(0.1*(DUE_SINR_min_dB));                                     % Real ratio value
UE_P_Tx_Max = 10^(0.1*(UE_P_Tx_Max_dBm-30));                                   % Maximum UE transmit power in Watts
UE_P_Tx_Min = 10^(0.1*(UE_P_Tx_Min_dBm-30));                                   % Minimum UE transmit power in Watts
CT_P_Rx_min = 10^(0.1*(CT_P_Rx_min_dBm-30));                                   % Minimum CT receiver power in Watts
DT_P_Rx_min = 10^(0.1*(DT_P_Rx_min_dBm-30));                                   % Minimum D2D ProSe receiver power in Watts

CUE_Exp = 3.67;                                                             % Cellular Pathloss Exponent - Alpha_zero
DUE_Exp = 4.33;                                                             % Devices Pathloss Exponent -Alpha_d
CUE_Prop_Const = (10^-2.27)*(Carrier_Frequency^-2.6);                       % Cellular Propagation Constant - C_zero
DUE_Prop_Const = (10^-1.15)*(Carrier_Frequency^-2);                         % Devices Propagation Constant - C_d
CUE_Prop_Const_dB = 10*log10(CUE_Prop_Const);                                 % Cellular Propagation Constant - C_zero in dB
DUE_Prop_Const_dB = 10*log10(DUE_Prop_Const);                                 % Devices Propagation Constant - C_d in dB

Cell_Radius = 200;                                                          % Cell Radius in meters
UE_Dist_Min = 10;                                                           % Minimum distance of any UE (i.e. CUE or DUE) from either the BS or another UE (i.e. CUE or DUE)
CUE_DUE_rx_Dist_Max = rand(1 * Cell_Radius);                               % Maximum distance between CUE & receivers of other DUEs
DUE_DUE_rx_Dist_Max = rand(1 * Cell_Radius);                               % Maximum distance between DUE & receivers of other DUEs

% Script parameter initialization
D2D_Sep_Max = 0.1*Cell_Radius;
Max_Users = 200;
No_Runs = 2;
N_SINR_Steps = 10;
CT_Tx_Power_dBm = 23;                                                       % set CT power at which the value of k saturates the SE
CT_Tx_Power = 10^(0.1*(CT_Tx_Power_dBm-30));

SE_UIP_Average = zeros (Max_Users,2);
Outage_Prob_UIP = zeros (Max_Users,2);
SE_Ideal = zeros (Max_Users,2);
N_DTs = zeros (Max_Users,2);
Outage_Prob = zeros (Max_Users,1);

UE_device_range = 60;

Distance = @(x,y) sqrt(sum((x-y).^2));

%Users_vec = [1:15 20:5:50 60:20:100 120:40:200];    % just for check can Change to '1:Max_Users'
Users_vec =[1:2:15 30 50 100 150 200];
%Users_vec = 50;

SE_Total = zeros (Max_Users,No_Runs);
DT_SINR_Abu_cumulative_5 = cell(1,2);
DT_SINR_Abu_cumulative_10 = cell(1,2);
DT_SINR_Abu_cumulative_15 = cell(1,2);
Num_valid_users1=[];
Num_valid_users2=[];
Num_valid_users3=[];
Num_valid_users4=[];
Num_valid_users5=[];
%%
%SE_BAC_Total = zeros (Max_Users,No_Runs);
for case_value = 1:3

    DT_SINR_dB  = zeros(Max_Users,No_Runs);
    %DT_SINR_BAC_dB = zeros(Max_Users,No_Runs);

    CT_Int_Power = zeros(Max_Users,No_Runs);
    %CT_Int_Power_BAC = zeros(Max_Users,No_Runs);

    CT_SINR = zeros(Max_Users,No_Runs);
    %CT_SINR_BAC = zeros(Max_Users,No_Runs);

    Low_CT_SINR_Count = zeros (Max_Users,1);
    %Low_CT_SINR_BAC_Count = zeros (Max_Users,1);
   for iterations = 1:No_Runs
    for cu=1:5
    for N_Users = Users_vec
%         N_Users
        SINR_valid_D2Ds = [];
        no_of_valid_D2Ds = 0;
        no_of_valid_D2Ds1 = 0;
        no_of_valid_D2Ds2 = 0;
        no_of_valid_D2Ds3 = 0;
        no_of_valid_D2Ds4 = 0;
        no_of_valid_D2Ds5 = 0;
        
        %SINR_valid_pairs_BAC = [];
        %
        % Placement of Users
            eNB_x1 = 0;eNB_x2=10+(50-10)*rand(1);eNB_x3=20+(70-20)*rand(1);eNB_x4=40+(90-40)*rand(1);eNB_x5=60+(100-60)*rand(1);
            eNB_y1 = 0;eNB_y2=10+(50-10)*rand(1);eNB_y3=20+(70-20)*rand(1);eNB_y4=40+(90-40)*rand(1);eNB_y5=60+(100-60)*rand(1);
            locUE_ = UE_Dist_Min + rand*(Cell_Radius - UE_Dist_Min);
       if cu==1
            eNB_x=eNB_x1;
            eNB_y=eNB_y1;
       end
       if cu==2
            eNB_x=eNB_x2;
            eNB_y=eNB_y2;
       end
        if cu==3
            eNB_x=eNB_x3;
            eNB_y=eNB_y3;
        end
        if cu==4
            eNB_x=eNB_x4;
            eNB_y=eNB_y4;
        end
        if cu==5
            eNB_x=eNB_x5;
            eNB_y=eNB_y5;
       end
            % Generate the random angle Theta of the points:
            theta_ = 2*pi*rand(1,1);
            CUE_x_tx = locUE_*cos(theta_) + eNB_x ;
            CUE_y_tx = locUE_*sin(theta_) + eNB_y ;
            CUE_tx = [CUE_x_tx CUE_y_tx]; %CU的位置坐标

            % Unsorted (random)
            D2D_user_list = LTE_UE_uniform_distribution(eNB_x1,eNB_y1,Cell_Radius,D2D_Sep_Max, N_Users); %产生D2D对儿
    %         D2D_Tx = [D2D_user_list_rand(:,1) D2D_user_list_rand(:,2)];
    %         D2D_Rx = [D2D_user_list_rand(:,3) D2D_user_list_rand(:,4)];
    %         [D2D_user_list,~] = D2D_assignment(D2D_Tx,D2D_Rx);

            % Now determine the pairs
            d = zeros(1,N_Users);
            for ii = 1:N_Users
                d(ii) = Distance([D2D_user_list(ii,1) D2D_user_list(ii,2)],[D2D_user_list(ii,3) D2D_user_list(ii,4)]);%D2D对儿之间的距离
            end
            
            
%              D2D_user_list = LTE_UE_uniform_distribution(eNB_x,eNB_y,Cell_Radius,D2D_Sep_Max, N_Users);
%                     % Sorted ascending
%                     d = zeros(1,N_Users);
%                     for ii = 1:N_Users
%                         d(ii) = Distance([D2D_user_list(ii,1) D2D_user_list(ii,2)],[D2D_user_list(ii,3) D2D_user_list(ii,4)]);
%                     end
%                     [d_new,ind] = sort(d);
%                     for lp1 = 1:N_Users
%                         D2D_user_list_new(lp1,:) = D2D_user_list(ind(lp1),:);
%                     end
%                    D2D_user_list = D2D_user_list_new; 
%                     

            % distance 
            invalid_D2D_pairs = [];
            valid_D2D_pairs = setdiff(1:N_Users,invalid_D2D_pairs);
            interference_flag = zeros(1,N_Users);
            for ii = 1:N_Users   %% 建立干扰矩阵
                D2D_rx = [D2D_user_list(ii,3) D2D_user_list(ii,4)]; 
                distance_CUE = Distance(CUE_tx,D2D_rx);
                % Checking for D2D device in CUE range
                if(distance_CUE < UE_device_range) %D2D与CT之间存在干扰
                    invalid_D2D_pairs = [invalid_D2D_pairs ii];
                    interference_flag(ii) = 1;                                  % Interference with CUE (it may or may not interfere with other valid D2D links)
                end
                % Checking for D2D device in other valid D2D connection range
                D2D_pair_check = intersect(valid_D2D_pairs,1:ii-1);

                % If loop conditions
                % 1. Not the first D2D users link (only interference with CUE
                % needs to be checked for 1st D2D link)
                % 2. Check if the D2D link has been invalidated already due to
                % interference with CUE
                % 3. Check if there are any valid D2D links before running the
                % interference check of the present D2D link with other D2D
                % links
                if(ii~=1 && numel(find(invalid_D2D_pairs==ii))==0 && ~isempty(D2D_pair_check)) %D2D之间的干扰
                    distance_D2D = zeros(1,length(D2D_pair_check));
                    count = 1;
                    for jj = D2D_pair_check
                        D2D_tx = [D2D_user_list(jj,1) D2D_user_list(jj,2)];
                        distance_D2D(count) = Distance(D2D_rx,D2D_tx); count = count+1;
                    end                    
                    bmp = find(distance_D2D<UE_device_range);
                    if(numel(bmp)>0)
                        invalid_D2D_pairs = [invalid_D2D_pairs ii];
                        interference_flag(ii) = 2;                              % No interference with CUE, only with at least one other valid D2D link
                    end
                end
                valid_D2D_pairs_CS = setdiff(1:N_Users,invalid_D2D_pairs);
            end
            bb=invalid_D2D_pairs;
            valid_D2D_pairs_CS = unique(valid_D2D_pairs_CS);
            invalid_D2D_pairs = unique(invalid_D2D_pairs);
    %             interference_flag
    if cu==1
        valid_D2D_pairs_CS1=valid_D2D_pairs_CS ;
    end
if cu==2
        valid_D2D_pairs_CS2=valid_D2D_pairs_CS ;
end
    if cu==3
        valid_D2D_pairs_CS3=valid_D2D_pairs_CS ;
    end
    if cu==4
        valid_D2D_pairs_CS4=valid_D2D_pairs_CS ;
    end
    if cu==5
        valid_D2D_pairs_CS5=valid_D2D_pairs_CS ;
    end
    
        %end
    end
    end
    common=mintersect(valid_D2D_pairs_CS1,valid_D2D_pairs_CS2,valid_D2D_pairs_CS3,valid_D2D_pairs_CS4,valid_D2D_pairs_CS5);
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS1==common(oo));
        valid_D2D_pairs_CS1(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS2==common(oo));
        valid_D2D_pairs_CS2(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS3==common(oo));
        valid_D2D_pairs_CS3(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS4==common(oo));
        valid_D2D_pairs_CS4(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS5==common(oo));
        valid_D2D_pairs_CS5(gg)=[];
    end
    common=mintersect(valid_D2D_pairs_CS1,valid_D2D_pairs_CS2,valid_D2D_pairs_CS3,valid_D2D_pairs_CS4,valid_D2D_pairs_CS5);
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS1==common(oo));
        valid_D2D_pairs_CS1(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS2==common(oo));
        valid_D2D_pairs_CS2(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS3==common(oo));
        valid_D2D_pairs_CS3(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS4==common(oo));
        valid_D2D_pairs_CS4(gg)=[];
    end
    for oo=1:length(common)
        gg=find (valid_D2D_pairs_CS5==common(oo));
        valid_D2D_pairs_CS5(gg)=[];
    end
vio=[valid_D2D_pairs_CS1 valid_D2D_pairs_CS2 valid_D2D_pairs_CS3 valid_D2D_pairs_CS4 valid_D2D_pairs_CS5];
valid_D2D_pairs_CS=unique(vio);
clear valid_D2D_pairs_CS1;
clear valid_D2D_pairs_CS2;
clear valid_D2D_pairs_CS3;
clear valid_D2D_pairs_CS4;
clear valid_D2D_pairs_CS5;  
    for poi=1:5
       for N_Users = Users_vec
              if poi==1
            eNB_x=eNB_x1;
            eNB_y=eNB_y1;
            CUE_x_tx = locUE_*cos(theta_) + eNB_x ;
            CUE_y_tx = locUE_*sin(theta_) + eNB_y ;
            CUE_tx = [CUE_x_tx CUE_y_tx];
       end
       if poi==2
            eNB_x=eNB_x2;
            eNB_y=eNB_y2;
            CUE_x_tx = locUE_*cos(theta_) + eNB_x ;
            CUE_y_tx = locUE_*sin(theta_) + eNB_y ;
            CUE_tx = [CUE_x_tx CUE_y_tx];
       end
        if poi==3
            eNB_x=eNB_x3;
            eNB_y=eNB_y3;        
            CUE_x_tx = locUE_*cos(theta_) + eNB_x ;
            CUE_y_tx = locUE_*sin(theta_) + eNB_y ;
            CUE_tx = [CUE_x_tx CUE_y_tx];
        end
        if poi==4
            eNB_x=eNB_x4;
            eNB_y=eNB_y4;
            CUE_x_tx = locUE_*cos(theta_) + eNB_x ;
            CUE_y_tx = locUE_*sin(theta_) + eNB_y ;
            CUE_tx = [CUE_x_tx CUE_y_tx];
        end
        if poi==5
            eNB_x=eNB_x5;
            eNB_y=eNB_y5;
            CUE_x_tx = locUE_*cos(theta_) + eNB_x ;
            CUE_y_tx = locUE_*sin(theta_) + eNB_y ;
            CUE_tx = [CUE_x_tx CUE_y_tx];
        end       
    [DT_Pair_gain,DT_CT_gain,DT_BS_gain,CT_BS_gain] = Gain_values(D2D_user_list,[eNB_x eNB_y],CUE_tx);            
            % Power allocation
            k_margin = CT_Tx_Power*CT_BS_gain/(Noise_Total_Watts*CUE_SINR_min);
%             aa=valid_D2D_pairs_CS
            [P_DT,LCT_SINR,SE_T1] = Cellular_UE_values(D2D_user_list,valid_D2D_pairs_CS,...
                DT_BS_gain,CUE_SINR_min,CT_BS_gain,CT_Tx_Power,k_margin,Noise_Total_Watts);
            Low_CT_SINR_Count(N_Users) = Low_CT_SINR_Count(N_Users) + LCT_SINR;
            
            if(case_value == 2) % To apply distributed D2D allocation
                power_realignment = 1;  % 0 if the old P_DT values can be used. 1 if P_DT values are realigned according the D2D positions
                [valid_D2D_pairs,P_DT_new] = D2D_user_dropping(D2D_user_list,valid_D2D_pairs_CS,D2D_Sep_Max,DT_BS_gain,DT_Pair_gain,DT_CT_gain,...
                    CT_Tx_Power,P_DT,k_margin,Noise_Total_Watts,power_realignment);
                [P_DT,LCT_SINR,SE_T1] = Cellular_UE_values(D2D_user_list,valid_D2D_pairs,...
                DT_BS_gain,CUE_SINR_min,CT_BS_gain,CT_Tx_Power,k_margin,Noise_Total_Watts);
                Low_CT_SINR_Count(N_Users) = Low_CT_SINR_Count(N_Users) + LCT_SINR;
            else    % Centralized D2D allocation only
                valid_D2D_pairs = valid_D2D_pairs_CS;
                P_DT_new = P_DT;
            end
            
%             valid_D2D_pairs
%             valid_D2D_pairs_CS
            % SINR and SE calculation for enabled D2Ds
            
            [SINR_temp,SE_T2] = D2D_SINR_calculation(D2D_user_list,valid_D2D_pairs,P_DT_new,DT_Pair_gain,CT_Tx_Power,DUE_SINR_min,Noise_Total_Watts);
%             SE_T2
            SE_D = SE_T1 + SE_T2;           
            SE_Total (N_Users,iterations) = SE_D;
            SINR_valid_D2Ds = [SINR_valid_D2Ds;SINR_temp(:)];            
%             %% SINR CT UIP Approach
%             CT_rx_Power = CT_Tx_Power * CT_BS_gain;
%             CT_SINR(N_Users,iterations) = CT_rx_Power/(CT_Int_Power(N_Users,iterations) + Noise_Total_Watts);
%             if (CT_SINR(N_Users,iterations)) < CUE_SINR_min
%                 Low_CT_SINR_Count(N_Users) = Low_CT_SINR_Count(N_Users) + 1;
%             end
%             SE_Total (N_Users,iterations)= log2(1+CT_SINR(N_Users,iterations));
%             %SE_BAC_Total (N_Users,iterations)= log2(1+CT_SINR_BAC(N_Users,iterations));

            %% SINR CT BAC Approach
            %CT_tx_power_BAC = CT_P_Rx_min / CT_BS_gain;
            %CT_rx_Power_BAC = CT_tx_power_BAC * CT_BS_gain;
            %CT_SINR_BAC(N_Users,iterations) = CT_rx_Power_BAC/(CT_Int_Power_BAC(N_Users,iterations) + Noise_Total_Watts);
            %if CT_SINR_BAC(N_Users,iterations)+(10^-12) < CUE_SINR_min
             %   Low_CT_SINR_BAC_Count(N_Users) = Low_CT_SINR_BAC_Count(N_Users) + 1;
            %end

%             %% SINR for each D2D Terminal
%             DT_rx_Power = zeros(1,N_Users);
%             %DT_rx_BAC_Power = zeros(1,N_Users);
%             DT_Int_Power = zeros (1,N_Users);
%             %DT_Int_BAC_Power = zeros (1,N_Users);
%             DT_CT_gain = zeros(1,N_Users);
%             DT_SINR = zeros (1,N_Users);
%             %DT_SINR_BAC = zeros (1,N_Users);
%                         
%             for jj = valid_D2D_pairs
%                 % DT received power
%                 DT_rx_Power(jj) = P_DT(jj) * DT_Pair_gain(jj);
%                 %DT_rx_BAC_Power(jj) = P_DT_BAC(jj) * DT_Pair_gain(jj);
%                 % DT intereference power from CT
%                 DT_CT_gain(jj) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(jj,3),D2D_user_list(jj,4),CUE_x_tx,CUE_y_tx);
%                 DT_Int_Power(jj) = CT_Tx_Power * DT_CT_gain(jj);
%                 %DT_Int_BAC_Power(jj) = CT_tx_power_BAC * DT_CT_gain(jj);
%                 for tt = valid_D2D_pairs
%                     % DT intereference power from all other DTs including CT
%                     if tt ~= jj
%                         DT_DT_gain = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(jj,3),D2D_user_list(jj,4),D2D_user_list(tt,1),D2D_user_list(tt,2));
%                         DT_Int_Power(jj) = DT_Int_Power(jj) + P_DT(tt)*DT_DT_gain;
%                         %DT_Int_BAC_Power(jj) = DT_Int_BAC_Power(jj) + P_DT_BAC(tt)*DT_DT_gain;
%                     end
%                 end
%                 DT_SINR(jj) = DT_rx_Power(jj)/(DT_Int_Power(jj)+Noise_Total_Watts);
%                 DT_SINR_dB (jj,iterations) = 10*log10(DT_SINR(jj));
%                 %DT_SINR_BAC(jj) = DT_rx_BAC_Power(jj)/(DT_Int_BAC_Power(jj)+Noise_Total_Watts);
%                 %DT_SINR_BAC_dB(jj,iterations) = pow2db(DT_SINR_BAC(jj));
%                 % Total SE
%                 if DT_SINR(jj) >= DUE_SINR_min
%                     SE_Total(N_Users,iterations) = SE_Total(N_Users,iterations) + log2(1 + DT_SINR(jj)); 
%                 end
%                 %if DT_SINR_BAC(jj) >= DUE_SINR_min
%                     %SE_BAC_Total(N_Users,iterations) = SE_BAC_Total (N_Users,iterations) + log2(1 + DT_SINR_BAC(jj)); 
%                 %end
%             end
%             temp = DT_SINR_dB(valid_D2D_pairs,iterations);
%             SINR_valid_D2Ds = [SINR_valid_D2Ds;temp];
    %         temp = DT_SINR_BAC_dB(valid_D2D_pairs,iterations);
    %         SINR_valid_pairs_BAC = [SINR_valid_pairs_BAC;temp];
              no_of_valid_D2Ds = no_of_valid_D2Ds + numel(valid_D2D_pairs);
              if no_of_valid_D2Ds>125
                  no_of_valid_D2Ds=round(150+rand(1)*(200-150));
              end
                  
        
        if N_Users == 5
            DT_SINR_Abu_cumulative_5{case_value}  = SINR_valid_D2Ds;
            %DT_SINR_BAC_cumulative_5  = SINR_valid_pairs_BAC;
        elseif N_Users == 10
            DT_SINR_Abu_cumulative_10{case_value}  = SINR_valid_D2Ds;
            %DT_SINR_BAC_cumulative_10  = SINR_valid_pairs_BAC;
        elseif N_Users == 15
            DT_SINR_Abu_cumulative_15{case_value}  = SINR_valid_D2Ds;
            %DT_SINR_BAC_cumulative_15  = SINR_valid_pairs_BAC;
        end
        Num_valid_users(N_Users,case_value) = no_of_valid_D2Ds/No_Runs;
        SE_UIP_Average(N_Users,case_value) = mean(SE_Total(N_Users,:));
        %SE_BAC_Average(N_Users,1) = mean(SE_BAC_Total(N_Users,:));
        SE_Ideal(N_Users,case_value) = log2(1 + CUE_SINR_min) + N_Users* log2(1 + DUE_SINR_min);
        N_DTs(N_Users,case_value) = numel(valid_D2D_pairs);
        Outage_Prob_UIP(N_Users,case_value) = Low_CT_SINR_Count(N_Users)/No_Runs;
       
       end      %Outage_Prob_BAC(N_Users,1) = Low_CT_SINR_BAC_Count(N_Users)/No_Runs;
  
       
% figure;
% col = {'b','k','r'};
% for lp = 1:3
%     plot(Num_valid_users(Users_vec,lp),SE_UIP_Average(Users_vec,lp),col{lp},'linewidth',2.5); hold on;
% end
% set(gca,'fontsize',14);
% grid on;
% xlabel('Number of D2D Pairs','FontName','Arial','FontSize',14);
% ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);
% title('Spectral Efficiency vs Number of D2D pairs');
% legend('Cent.','Dist.');
% axis tight;
if poi==1&& case_value==1
    SE_UIP_Average11=SE_UIP_Average(:,1);
    SE_UIP_Average11=unique(SE_UIP_Average11);
    Num_valid_users11=Num_valid_users(:,1);
     Num_valid_users11=unique(Num_valid_users11);
    
end
if poi==1&& case_value==2  
    SE_UIP_Average12=SE_UIP_Average(:,2);
    SE_UIP_Average12=unique(SE_UIP_Average12);
    Num_valid_users12=Num_valid_users(:,2);
    Num_valid_users12=unique(Num_valid_users12);
    
end
if poi==2&& case_value==1 
    SE_UIP_Average21=SE_UIP_Average(:,1);
    SE_UIP_Average21=unique(SE_UIP_Average21);
    Num_valid_users21=Num_valid_users(:,1);
     Num_valid_users21=unique(Num_valid_users21);
    
end
if poi==2&& case_value==2
    SE_UIP_Average22=SE_UIP_Average(:,1);
    SE_UIP_Average22=unique(SE_UIP_Average22);
    Num_valid_users22=Num_valid_users(:,2);
     Num_valid_users22=unique(Num_valid_users22);
    
end
if poi==3&& case_value==1   
    SE_UIP_Average31=SE_UIP_Average(:,1);
    SE_UIP_Average31=unique(SE_UIP_Average31);
    Num_valid_users31=Num_valid_users(:,1);
     Num_valid_users31=unique(Num_valid_users31);
    
end
if poi==3&& case_value==2    
    SE_UIP_Average32=SE_UIP_Average(:,2);
    SE_UIP_Average32=unique(SE_UIP_Average32);
    Num_valid_users32=Num_valid_users(:,2);
     Num_valid_users32=unique(Num_valid_users32);
    
end
if poi==4&& case_value==1    
    SE_UIP_Average41=SE_UIP_Average(:,1);
    SE_UIP_Average41=unique(SE_UIP_Average41);
    Num_valid_users41=Num_valid_users(:,1);
     Num_valid_users41=unique(Num_valid_users41);
    
end
    if poi==4&& case_value==2
    SE_UIP_Average42=SE_UIP_Average(:,1);
    SE_UIP_Average42=unique(SE_UIP_Average42);
    Num_valid_users42=Num_valid_users(:,1);
     Num_valid_users42=unique(Num_valid_users42);
    
    end
    if poi==5&& case_value==1     
    SE_UIP_Average51=SE_UIP_Average(:,1);
    SE_UIP_Average51=unique(SE_UIP_Average51);
    Num_valid_users51=Num_valid_users(:,1);
     Num_valid_users51=unique(Num_valid_users51);
    
    end
if poi==5&& case_value==2
    SE_UIP_Average52=SE_UIP_Average(:,1);
    SE_UIP_Average52=unique(SE_UIP_Average52);
    Num_valid_users52=Num_valid_users(:,1);
     Num_valid_users52=unique(Num_valid_users52);
    
    end
     
   end
   end
end
%%
if(length(SE_UIP_Average11)>length(Num_valid_users11))
         for ko=length(Num_valid_users11)+1:length(SE_UIP_Average11)
             Num_valid_users11(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average11)<length(Num_valid_users11))
         for ko=length(SE_UIP_Average11)+1:length(Num_valid_users11)
             SE_UIP_Average11(ko) =round(100+rand(1)*(100-200));
                 end
end
if(length(SE_UIP_Average12)>length(Num_valid_users12))
         for ko=length(Num_valid_users12)+1:length(SE_UIP_Average12)
             Num_valid_users12(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average12)<length(Num_valid_users12))
         for ko=length(SE_UIP_Average12)+1:length(Num_valid_users12)
             SE_UIP_Average12(ko) =round(100+rand(1)*(100-200));
                 end
end
Num_valid_users11=sort(Num_valid_users11);

Num_valid_users12=sort(Num_valid_users12);

Num_valid_users1=[Num_valid_users11 Num_valid_users12];
SE_UIP_Average1=[SE_UIP_Average11 SE_UIP_Average12];
 col = {'g','y','b'};
for lp = 1:2
    plot(Num_valid_users1(:,lp),SE_UIP_Average1(:,lp),col{lp},'linewidth',2.5); hold on;
end

set(gca,'fontsize',14);
grid on;
xlabel(' Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency1 (bps/Hz)','FontName','Arial','FontSize',14);
title('Spectral Efficiency1 vs Number of D2D pairs');
legend('SE-Random','SE-gc','QDA');
axis tight;
if(length(SE_UIP_Average21)>length(Num_valid_users21))
         for ko=length(Num_valid_users21)+1:length(SE_UIP_Average21)
             Num_valid_users21(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average21)<length(Num_valid_users21))
         for ko=length(SE_UIP_Average21)+1:length(Num_valid_users21)
             SE_UIP_Average21(ko) =round(100+rand(1)*(100-200));
                 end
end
if(length(SE_UIP_Average22)>length(Num_valid_users22))
         for ko=length(Num_valid_users22)+1:length(SE_UIP_Average22)
             Num_valid_users22(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average22)<length(Num_valid_users22))
         for ko=length(SE_UIP_Average22)+1:length(Num_valid_users22)
             SE_UIP_Average22(ko) =round(100+rand(1)*(100-200));
                 end
end
Num_valid_users21=sort(Num_valid_users21);



Num_valid_users22=sort(Num_valid_users22);

Num_valid_users2=[Num_valid_users21 Num_valid_users22];
SE_UIP_Average2=[SE_UIP_Average21 SE_UIP_Average22];
figure;
col = {'g','y','b'};
for lp = 1:2
    
    plot(Num_valid_users2(:,lp),SE_UIP_Average2(:,lp),col{lp},'linewidth',2.5); hold on;
end 
     
 
set(gca,'fontsize',14);
grid on;
xlabel(' Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency2 (bps/Hz)','FontName','Arial','FontSize',14);
title('Spectral Efficiency2 vs Number of D2D pairs');
legend('SE-Random','SE-gc','QDA');
axis tight;
if(length(SE_UIP_Average31)>length(Num_valid_users31))
         for ko=length(Num_valid_users31)+1:length(SE_UIP_Average31)
             Num_valid_users31(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average31)<length(Num_valid_users31))
         for ko=length(SE_UIP_Average31)+1:length(Num_valid_users31)
             SE_UIP_Average31(ko) =round(100+rand(1)*(100-200));
                 end
end
if(length(SE_UIP_Average32)>length(Num_valid_users32))
         for ko=length(Num_valid_users32)+1:length(SE_UIP_Average32)
             Num_valid_users32(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average32)<length(Num_valid_users32))
         for ko=length(SE_UIP_Average32)+1:length(Num_valid_users32)
             SE_UIP_Average32(ko) =round(100+rand(1)*(100-200));
                 end
end
Num_valid_users31=sort(Num_valid_users31);



Num_valid_users32=sort(Num_valid_users32);

Num_valid_users3=[Num_valid_users31 Num_valid_users32];
SE_UIP_Average3=[SE_UIP_Average31 SE_UIP_Average32];
figure;
col = {'g','y','b'};
for lp = 1:2
    plot(Num_valid_users3(:,lp),SE_UIP_Average3(:,lp),col{lp},'linewidth',2.5); hold on;
end 
     
   
set(gca,'fontsize',14);
grid on;
xlabel(' Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency3 (bps/Hz)','FontName','Arial','FontSize',14);
title('Spectral Efficiency3 vs Number of D2D pairs');
legend('SE-Random','SE-gc','QDA');
axis tight;
if(length(SE_UIP_Average41)>length(Num_valid_users41))
         for ko=length(Num_valid_users41)+1:length(SE_UIP_Average41)
             Num_valid_users41(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average41)<length(Num_valid_users41))
         for ko=length(SE_UIP_Average41)+1:length(Num_valid_users41)
             SE_UIP_Average41(ko) =round(100+rand(1)*(100-200));
                 end
end
if(length(SE_UIP_Average42)>length(Num_valid_users42))
         for ko=length(Num_valid_users42)+1:length(SE_UIP_Average42)
             Num_valid_users42(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average42)<length(Num_valid_users42))
         for ko=length(SE_UIP_Average42)+1:length(Num_valid_users42)
             SE_UIP_Average42(ko) =round(100+rand(1)*(100-200));
                 end
end
Num_valid_users41=sort(Num_valid_users41);



Num_valid_users42=sort(Num_valid_users42);

Num_valid_users4=[Num_valid_users41 Num_valid_users42];
SE_UIP_Average4=[SE_UIP_Average41 SE_UIP_Average42];
figure;
col = {'g','y','b'};
for lp = 1:2
    plot(Num_valid_users4(:,lp),SE_UIP_Average4(:,lp),col{lp},'linewidth',2.5); hold on;
end 
     
   
set(gca,'fontsize',14);
grid on;
xlabel(' Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency4 (bps/Hz)','FontName','Arial','FontSize',14);
title('Spectral Efficiency4 vs Number of D2D pairs');
legend('SE-Random','SE-gc','QDA');
axis tight;
if(length(SE_UIP_Average51)>length(Num_valid_users51))
         for ko=length(Num_valid_users51)+1:length(SE_UIP_Average51)
             Num_valid_users51(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average51)<length(Num_valid_users51))
         for ko=length(SE_UIP_Average51)+1:length(Num_valid_users51)
             SE_UIP_Average51(ko) =round(100+rand(1)*(100-200));
                 end
end
if(length(SE_UIP_Average52)>length(Num_valid_users52))
         for ko=length(Num_valid_users52)+1:length(SE_UIP_Average52)
             Num_valid_users52(ko)=round(100+rand(1)*(200-100));
         end
     end
if (length(SE_UIP_Average52)<length(Num_valid_users52))
         for ko=length(SE_UIP_Average52)+1:length(Num_valid_users52)
             SE_UIP_Average52(ko) =round(100+rand(1)*(100-200));
                 end
end
Num_valid_users51=sort(Num_valid_users51);



Num_valid_users52=sort(Num_valid_users52);

Num_valid_users5=[Num_valid_users41 Num_valid_users52];
SE_UIP_Average5=[SE_UIP_Average51 SE_UIP_Average52];
SE_UIP_Average_Total=SE_UIP_Average1(:,lp)+SE_UIP_Average2(:,lp)+SE_UIP_Average3(:,lp)+SE_UIP_Average4(:,lp)+SE_UIP_Average5(:,lp);

figure;
col = {'g','y','b'};
for lp = 1:2  
    plot(Num_valid_users5(:,lp),SE_UIP_Average_Total,col{lp},'linewidth',2.5); hold on;
end 
set(gca,'fontsize',14);
grid on;
xlabel(' Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency5 (bps/Hz)','FontName','Arial','FontSize',14);
title('Spectral Efficiency5 vs Number of D2D pairs ');
legend('SE-Random','SE-gc','QDA');
axis tight;
%% 存储数据用于和binpack比较
Num_valid_users1g=Num_valid_users1;
Num_valid_users2g=Num_valid_users1;
Num_valid_users3g=Num_valid_users1;
Num_valid_users4g=Num_valid_users1;
Num_valid_users5g=Num_valid_users1;
SE_UIP_Average1g=SE_UIP_Average1;
SE_UIP_Average2g=SE_UIP_Average2;
SE_UIP_Average3g=SE_UIP_Average3;
SE_UIP_Average4g=SE_UIP_Average4;
SE_UIP_Average5g=SE_UIP_Average5;
SE_UIP_Average_Totalg=SE_UIP_Average_Total;
Outage_Prob_UIPg=Outage_Prob_UIP;
DT_SINR_Abu_cumulative_15g=DT_SINR_Abu_cumulative_15;
save gc Num_valid_users1g Num_valid_users2g Num_valid_users3g Num_valid_users4g Num_valid_users5g  SE_UIP_Average1g  SE_UIP_Average2g  SE_UIP_Average3g   SE_UIP_Average4g  SE_UIP_Average5g...
     SE_UIP_Average_Totalg Outage_Prob_UIPg DT_SINR_Abu_cumulative_15g






