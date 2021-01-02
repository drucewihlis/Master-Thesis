% This function uniformly distributes n UEs in a hexagonal cell of radius r
function [UExy] = LTE_UE_uniform_distribution(center_x,center_y,areaRad,D2DRad, N_Users, sector)
% Defining BS (Cell Center) Location
% cell center
%global LTE_config
x_ = center_x;
y_ = center_y;
UE_Dist_Min = 10;   
% Radius of cell, i.e. sector 
r_=1*areaRad;
    if nargin < 6
        sector = 0;
    
    end
    UExy = zeros(N_Users,4);
    
    for i=1:N_Users
switch sector
    case 1
        
        % Cell/Hexagon vertices
        t_=linspace(0,2*pi,7);
        hexagonVertix_x = x_ + r_ * cos(t_);
        hexagonVertix_y = y_ + r_ * sin(t_); 
        
        % condition value to break while 
         label_='true';
         while strcmp(label_,'true')
             %generate the random radial distance point uniformly distributed between 0 and R
             locUE_ = UE_Dist_Min + (r_ - UE_Dist_Min)*sqrt(rand(1,1));
             % Generate the random angle Theta of the points:
             theta_=2*pi*rand(1,1);
             UE_x_ = locUE_*cos(theta_) + x_ ;
             UE_y_ = locUE_*sin(theta_) + y_ ;
          
             if inpolygon(UE_x_,UE_y_, hexagonVertix_x,hexagonVertix_y) 
                UE_x_tx =UE_x_ ;
                UE_y_tx =UE_y_ ;
                label_='false';
                break;
             end
             D2DRx_ = UE_Dist_Min + (D2DRad - UE_Dist_Min)*sqrt(rand(1,1));
            % Generate the random angle Theta of the points:
            theta_= 2*pi*rand(1,1);
            UE_x_rx = D2DRx_*cos(theta_) + UE_x_tx ;
            UE_y_rx = D2DRx_*sin(theta_) + UE_y_tx ;
            
            UExy (i,:) = [UE_x_tx,UE_y_tx,UE_x_rx,UE_y_rx];
             
          end
    case 2
        
         % Cell/Hexagon vertices
         hexagonVertix_x = x_ + r_ * cos((2:4)*pi/3);
         hexagonVertix_y = y_ + r_ * sin((2:4)*pi/3); 
         % condition value to break while 
         label_='true';
         while strcmp(label_,'true')
             %generate the random radial distance point uniformly distributed between 0 and R
              locUE_ = UE_Dist_Min + (r_ - UE_Dist_Min)*sqrt(rand(1,1));
              % Generate the random angle Theta of the points:
              theta_=2*pi*rand(1,1);
              UE_x_ = locUE_*cos(theta_) + x_ ;
              UE_y_ = locUE_*sin(theta_) + y_ ;
          
              if inpolygon(UE_x_,UE_y_, hexagonVertix_x,hexagonVertix_y) 
                  UE_x_tx =UE_x_ ;
                  UE_y_tx =UE_y_ ;
                  label_='false';
                  break;
              end
              D2DRx_ = UE_Dist_Min + (D2DRad - UE_Dist_Min)*sqrt(rand(1,1));
            % Generate the random angle Theta of the points:
            theta_= 2*pi*rand(1,1);
            UE_x_rx = D2DRx_*cos(theta_) + UE_x_tx ;
            UE_y_rx = D2DRx_*sin(theta_) + UE_y_tx ;
            
            UExy (i,:) = [UE_x_tx,UE_y_tx,UE_x_rx,UE_y_rx];
          end
  
    
    
    case 3
        % Cell/Hexagon vertices
        hexagonVertix_x = x_ + r_ * cos((4:6)*pi/3);
        hexagonVertix_y = y_ + r_ * sin((4:6)*pi/3); 
        
        % condition value to break while 
        label_='true';
        while strcmp(label_,'true')
            %generate the random radial distance point uniformly distributed between 0 and R
            locUE_ = UE_Dist_Min + (r_ - UE_Dist_Min)*sqrt(rand(1,1));
            % Generate the random angle Theta of the points:
            theta_=2*pi*rand(1,1);
            UE_x_ = locUE_*cos(theta_) + x_ ;
            UE_y_ = locUE_*sin(theta_) + y_ ;
            if inpolygon(UE_x_,UE_y_, hexagonVertix_x,hexagonVertix_y) 
              UE_x_tx =UE_x_ ;
              UE_y_tx =UE_y_ ;
              label_='false';
              break;
            end
            D2DRx_ = UE_Dist_Min + (D2DRad - UE_Dist_Min)*sqrt(rand(1,1));
            % Generate the random angle Theta of the points:
            theta_= 2*pi*rand(1,1);
            UE_x_rx = D2DRx_*cos(theta_) + UE_x_tx ;
            UE_y_rx = D2DRx_*sin(theta_) + UE_y_tx ;
            
            UExy (i,:) = [UE_x_tx,UE_y_tx,UE_x_rx,UE_y_rx];  
        end
    otherwise
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


