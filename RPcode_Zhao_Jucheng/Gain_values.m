function [DT_DT_gain,DT_CT_gain,DT_BS_gain,CT_BS_gain] = Gain_values(D2D_user_list,eNB_pos,CUE_pos)

N_D2D_pairs = size(D2D_user_list,1);
DT_CT_gain = zeros(1,N_D2D_pairs);
DT_BS_gain = zeros(1,N_D2D_pairs);
DT_DT_gain = zeros(N_D2D_pairs);

% Cellular UE position
CUE_x = CUE_pos(1);
CUE_y = CUE_pos(2);

% eNB position
eNB_x = eNB_pos(1);
eNB_y = eNB_pos(2);

for ii = 1:N_D2D_pairs
    for jj = 1:1:N_D2D_pairs
        % Gain from the j-th transmitter to the i-th receiver among the D2D pairs 
        DT_DT_gain(ii,jj) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(ii,3),D2D_user_list(ii,4),D2D_user_list(jj,1),D2D_user_list(jj,2));
    end
    
    % Gain from the i-th D2D receiver and the CT
    DT_CT_gain(1,ii) = LTE_channel_model_indoor_hotspot_NLOS(D2D_user_list(ii,3),D2D_user_list(ii,4),CUE_x,CUE_y);
    DT_BS_gain(1,ii) = LTE_channel_model_urban_micro_NLOS(D2D_user_list(ii,1),D2D_user_list(ii,2),eNB_x,eNB_y);
end

CT_BS_gain = LTE_channel_model_urban_micro_NLOS(CUE_x,CUE_y,eNB_x,eNB_y);