function Gain = LTE_channel_model_indoor_hotspot_NLOS( source_x,source_y, dest_x,dest_y)
% this is channel model which apply to device to device communication
%   Function calculates the pathloss in dB between two DUEs located at (dUE_sx,dUE_sy) and (dUE_dx,dUE_dy)
%   a carrier frequency fc is in GHz
%   The distribution of the shadow fading is log-normal with standard deviation sigma
%   Function follows the 3GPP TR 36.814 V9.0.0 Standard for Indoor Hotspot
%   NLOS applicable for 10m < d_D2D < 150m and Sigma = 4
%%
   %global LTE_config;
     %fc=LTE_config.frequency;
     %sigma=LTE_config.shadowing;
     fc=2.6;
     sigma=4; %shadow fading
     
%%     
   % Here distance between two devices is calcualted
      distance_ = pdist([source_x source_y;dest_x dest_y]);
   % one random number is generated for the log-normal fading of the channel 
      %rndnum=  normrnd(0,sigma);
      rndnum=  sigma*randn;
   %  total loss between two devices is calculated    
      loss = 43.3*log10(distance_) + 11.5 + 20*log10(fc) +rndnum ;
      Gain = 1/(db2pow(loss));
end
