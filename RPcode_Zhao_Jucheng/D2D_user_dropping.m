function [valid_D2D_pairs,P_DT] = D2D_user_dropping(D2D_user_list,valid_D2D_pairs,D2D_range,DT_BS_gain,DT_DT_gain,DT_CT_gain,...
    CT_Tx_Power,P_DT,k_margin,noise_power,power_realignment)

Distance = @(x,y) sqrt(sum((x-y).^2));  % Calculate distance between two positions

N_D2D_pairs = size(D2D_user_list,1);
Distance_Rx = zeros(N_D2D_pairs);
Distance_Tx = zeros(N_D2D_pairs);

for ii = 1:N_D2D_pairs
    for jj = 1:1:N_D2D_pairs
        Distance_Rx(ii,jj) = Distance(D2D_user_list(ii,[1 2]),D2D_user_list(jj,[3 4])); % Distance between i-th transmitter and j-th receiver
        Distance_Tx(ii,jj) = Distance(D2D_user_list(ii,[1 2]),D2D_user_list(jj,[1 2])); % Distance between i-th transmitter and j-th transmitter
    end
end

D2D_pairs_to_check = valid_D2D_pairs;   % The valid D2D pairs are the ones to be checked sequentially

% D2D_flags = zeros(1,N_D2D_pairs);
% replacement_D2D = zeros(1,N_D2D_pairs);

while(~isempty(D2D_pairs_to_check))
    k = D2D_pairs_to_check(1);                                  % Select the first D2D pair to check
    invalid_pairs = setdiff(1:N_D2D_pairs, valid_D2D_pairs);    % Find out all invalid pairs
    Rxs_in_range = find(Distance_Rx(k,:)<D2D_range);            % Find D2D Rxs in range
    Txs_in_range = find(Distance_Tx(k,:)<D2D_range);            % Find D2D Txs in range
    D2Ds_in_range = union(Rxs_in_range,Txs_in_range);           % A union of D2D Rxs and Txs in range gives the D2D pairs in range
    D2Ds_in_range = setdiff(D2Ds_in_range,k);                   % Remove the already associated D2D receiver
    D2Ds_in_range = intersect(D2Ds_in_range,invalid_pairs);     % Check for D2Ds that aren't activated within D2D range
    
    if(power_realignment == 1)
        P_DT = zeros(1,N_D2D_pairs);
        DT_BS_valid_pairs = DT_BS_gain(valid_D2D_pairs);
        P_DT_valid_pairs = ((k_margin-1)*noise_power/numel(valid_D2D_pairs))./DT_BS_valid_pairs;
        P_DT(valid_D2D_pairs) = P_DT_valid_pairs;
    else
%         P_DT_valid_pairs = P_DT(valid_D2D_pairs);
        % Do nothing
    end
    
    if(isempty(D2Ds_in_range))                                  % If there are no inactive D2Ds in range, do nothing and move on to the next UE
        D2D_pairs_to_check = D2D_pairs_to_check(2:end);
        valid_D2D_pairs = union(valid_D2D_pairs,k);
    else  % If there are other D2Ds in range that are not activated, check how much SE the non-activated ones offer. Choose the D2D pair that offers the highest SE
        
        rest_of_valid_pairs = setdiff(valid_D2D_pairs,k);
        DT_pair_self = DT_DT_gain(k,k);
        P_self = P_DT(k);
        DT_rx_power_self = P_self*DT_pair_self;
        int_D2D_self = DT_DT_gain(k,rest_of_valid_pairs);
        DT_int_power_self = sum(P_DT(rest_of_valid_pairs).*int_D2D_self) + CT_Tx_Power * DT_CT_gain(k);
        SINR_self = DT_rx_power_self/(DT_int_power_self + noise_power);
        SE_self = log2(1+SINR_self);    % Evaluate the SE of the D2D pair under consideration
        
        % Calculate of the SE of the D2D pairs in range
        DT_int_power = zeros(1,length(D2Ds_in_range));
        DT_rx_power = zeros(1,length(D2Ds_in_range));
        for lp = 1:length(D2Ds_in_range)
            int_D2D_gains = DT_DT_gain(D2Ds_in_range(lp),rest_of_valid_pairs);
            DT_int_power(lp) = sum(P_DT(rest_of_valid_pairs).*int_D2D_gains) + CT_Tx_Power * DT_CT_gain(D2Ds_in_range(lp));
            if(power_realignment == 1)
                new_valid_pairs = [D2Ds_in_range(lp) rest_of_valid_pairs];
                DT_BS_new_pairs = DT_BS_gain(new_valid_pairs);
                P_val = ((k_margin-1)*noise_power/numel(new_valid_pairs))./DT_BS_new_pairs;
                DT_rx_power(lp) = P_val(1)*DT_DT_gain(D2Ds_in_range(lp),D2Ds_in_range(lp));
            else
                DT_rx_power(lp) = P_DT(k);
            end
        end
        SINR_D2Ds_in_range = DT_rx_power./(DT_int_power + noise_power);
        SE_D2Ds_in_range = log2(1+SINR_D2Ds_in_range);  % Evaluate the SE of the D2Ds in range
        
        better_D2Ds = find(SE_D2Ds_in_range>SE_self);   % find the list of better D2Ds in terms of SE
        if(~isempty(better_D2Ds))   % If there are better D2Ds in terms of SE, replace the D2D with the best SE instead of the activated one in the valid D2Ds list
            [~,idx] = max(SE_D2Ds_in_range(better_D2Ds));
            D2D_pairs_to_check(1) = D2Ds_in_range(better_D2Ds(idx));
%             if(replacement_D2D(D2Ds_in_range(better_D2Ds(idx))) == k)
%                 D2D_pairs_to_check = D2D_pairs_to_check(2:end);
%                 continue;
%             end
%             replacement_D2D(k) = D2Ds_in_range(better_D2Ds(idx));
            if(power_realignment == 0)
                P_DT(D2Ds_in_range(better_D2Ds(idx))) = P_DT(k);
            end
            valid_D2D_pairs = setdiff(valid_D2D_pairs,k);
            valid_D2D_pairs = union(valid_D2D_pairs,D2Ds_in_range(better_D2Ds(idx)));
        else    % If none of the other D2Ds in range are better, don't make any changes and move on to the next activated D2D
            D2D_pairs_to_check = D2D_pairs_to_check(2:end);
            valid_D2D_pairs = union(valid_D2D_pairs,k);
        end
    end
end