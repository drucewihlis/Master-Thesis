function Gain = LTE_channel_model_urban_micro_NLOS(source_x,source_y, dest_x,dest_y)
%Cellular Channel Model
%   Function calculates the pathloss in dB between the CUEs and eNB located at (x_s,y_s) and (x_o,y_o) respectively and operating at a carrier frequency Fc (in GHz).
%   Underscore-s and underscore-o in the position coordinates denote the source and observation points respectively
%   The distribution of the shadow fading is log-normal with standard deviation Sigma
%   Function follows the 3GPP TR 36.814 V9.0.0 Standard for Urban Micro non-LOS applicable for 10m < d_CUE < 2,000m and Sigma = 3
%%
%      global LTE_config;
%      fc=LTE_config.frequency;
%      sigma=LTE_config.shadowing;
     fc=2.6;
     sigma=4; %shadow fading
%%
      % Here distance between two devices is calcualted
        distance_ = pdist([source_x source_y ; dest_x dest_y]);
       %disp(distance_);
        % one random number is generated for the log-normal fading of the channel 
        %rndnum= normrnd(0,sigma);
        rndnum=  sigma*randn;
       % total loss between two devices is calculated    
        loss = 36.7*log10(distance_) + 22.7 + 26*log10(fc) + rndnum;
        Gain = 1/(db2pow(loss));
end

