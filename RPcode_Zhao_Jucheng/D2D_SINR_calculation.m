function [SINR_valid_D2Ds,SE_Total] = D2D_SINR_calculation(D2D_user_list,valid_D2D_pairs,P_DT,DT_Pair_gain,CT_Tx_Power,DUE_SINR_min,Noise_Total_Watts)

N_Users = size(D2D_user_list,1);

DT_rx_Power = zeros(1,N_Users);
DT_Int_Power = zeros (1,N_Users);
DT_CT_gain = zeros(1,N_Users);
DT_SINR = zeros (1,N_Users);
SE_Total = 0;
DT_SINR_dB = zeros(1,N_Users);

for jj = valid_D2D_pairs
    % DT received power
    DT_rx_Power(jj) = P_DT(jj) * DT_Pair_gain(jj,jj);
    % DT intereference power from CT
    DT_Int_Power(jj) = CT_Tx_Power * DT_CT_gain(jj);
    %DT_Int_BAC_Power(jj) = CT_tx_power_BAC * DT_CT_gain(jj);
    for tt = valid_D2D_pairs
        % DT intereference power from all other DTs including CT
        if tt ~= jj
            DT_Int_Power(jj) = DT_Int_Power(jj) + P_DT(tt)*DT_Pair_gain(jj,tt);
        end
    end
    DT_SINR(jj) = DT_rx_Power(jj)/(DT_Int_Power(jj)+Noise_Total_Watts);
    DT_SINR_dB (jj) = 10*log10(DT_SINR(jj));
    % Total SE
    if DT_SINR(jj) >= DUE_SINR_min
        SE_Total = SE_Total + log2(1 + DT_SINR(jj));
    end
end
SINR_valid_D2Ds = DT_SINR_dB(valid_D2D_pairs);