% This function uniformly distributes n UEs in a hexagonal cell of radius r
function [UExy] = LTE_UE_uniform_distribution_upd(center_x,center_y,areaRad,D2DRad, N_Users, sector)
%%
% Defining BS (Cell Center) Location
% cell center
%global LTE_config
x_ = center_x;
y_ = center_y;
UE_Dist_Min = 10;   
% Radius of cell, i.e. sector 
r_=1*areaRad;
    if nargin < 6 %sector is not defined
       sector = 0; %case: otherwise
    end
    
    UExy = zeros(N_Users,4);
    
    for i=1:N_Users
        if sector==0
            %3 cases are deleted

            %generate the random radial distance point uniformly distributed between 0 and R
            locUE_ = UE_Dist_Min + (r_ - UE_Dist_Min)*sqrt(rand(1,1));
            % Generate the random angle Theta of the points:
            theta_= 2*pi*rand(1,1);
            UE_x_tx = locUE_*cos(theta_) + x_ ;
            UE_y_tx = locUE_*sin(theta_) + y_ ;
            
            D2DRx_ = UE_Dist_Min + (D2DRad - UE_Dist_Min)*sqrt(rand(1,1));
            % Generate the random angle Theta of the points:
            theta_= 2*pi*rand(1,1);
            UE_x_rx = D2DRx_*cos(theta_) + UE_x_tx ;
            UE_y_rx = D2DRx_*sin(theta_) + UE_y_tx ;
            
            UExy (i,:) = [UE_x_tx,UE_y_tx,UE_x_rx,UE_y_rx];           
        end       
    end   
end        


