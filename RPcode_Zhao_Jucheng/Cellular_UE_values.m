% Function that calculates the Tx power, interference power and SINR of
% cellular UE and thereby the SE from the cellular UE. It also evaluates if
% the cellular UE has gone into outage based on the SINR
function [P_DT,Low_CT_SINR_Count,SE_Total] = Cellular_UE_values(D2D_user_list,valid_D2D_pairs,...
    DT_BS_gain,CUE_SINR_min,CT_BS_gain,CT_Tx_Power,k_margin,Noise_Total_Watts)
CT_rx_Power = CT_Tx_Power * CT_BS_gain;
N_Users = size(D2D_user_list,1);
P_DT = zeros(1,N_Users);
CT_Int_Power = 0;
Low_CT_SINR_Count = 0;
%P_DT_BAC = zeros(1,N_Users);

for jj=valid_D2D_pairs                                                                   % User valid_D2D_pairs instead of 1:N_Users
    P_DT(jj) = ((k_margin-1)*Noise_Total_Watts/numel(valid_D2D_pairs))/DT_BS_gain(jj);           % Even share of the interference power to BS
    % DT_Tx power set by channel inversion based on the constant minimum DT receiver sensitivity
    %             P_DT_BAC(jj) = min(UE_P_Tx_Max, DT_P_Rx_min/DT_Pair_gain(jj));

    %Interference power for D2D pairs sharing the same RB as
    CT_Int_Power = CT_Int_Power + P_DT(jj)*DT_BS_gain(jj);
    %CT_Int_Power_BAC(N_Users,iterations) = CT_Int_Power_BAC(N_Users,iterations) + P_DT_BAC(jj)*DT_BS_gain(jj);
    CT_SINR = CT_rx_Power/(CT_Int_Power + Noise_Total_Watts);
if (CT_SINR) > CUE_SINR_min
    Low_CT_SINR_Count = 1;
    SE_Total = log2(1+CT_SINR);
else 
    break
end
end
%SE_BAC_Total (N_Users,iterations)= log2(1+CT_SINR_BAC(N_Users,iterations));