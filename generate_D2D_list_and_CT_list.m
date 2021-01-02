%saves 2 lists mapped to grid
clc;clear;close all;format compact;

N=30;
Cell_Radius = 110;
D2D_Sep_Max = 0.1*Cell_Radius;
Max_Users = 50;
CUE_SINR_min_dB	= 33; 
DUE_SINR_min_dB	= 46;
CUE_SINR_min = db2pow(CUE_SINR_min_dB); 
DUE_SINR_min = db2pow(DUE_SINR_min_dB); 
eNB1_x = -272; eNB1_y = -645;  %BSs in blue area
eNB2_x = -374; eNB2_y = -809;
eNB3_x = -491; eNB3_y = -662;
grid_Kazan=xlsread('Points.xlsx');
    eNB_img_x=-384; eNB_img_y = -697; 
    Cell_Radius_img= 210;
    Max_Users_single = 100;
    
D2D_user_list_all = LTE_UE_uniform_distribution_upd(eNB_img_x,eNB_img_y,Cell_Radius_img,D2D_Sep_Max, Max_Users);
[D2D_user_list_all_mapped,grid_Kazan_upd] = mapping_to_grid(grid_Kazan,D2D_user_list_all);
k=1;
for i=1:length(D2D_user_list_all_mapped)
    D2D_user_list_all_mapped_1col(k,:)=[D2D_user_list_all_mapped(i,1) D2D_user_list_all_mapped(i,2)];
    k=k+1;
    D2D_user_list_all_mapped_1col(k,:)=[D2D_user_list_all_mapped(i,3) D2D_user_list_all_mapped(i,4)];
    k=k+1;
end %D2D_user_list_all_mapped_1col -> OFDMPlanning to get Loss

[CUE1_x, CUE1_y] = generate_CT(eNB1_x,eNB1_y,Cell_Radius);
[CUE2_x, CUE2_y] = generate_CT(eNB2_x,eNB2_y,Cell_Radius);
[CUE3_x, CUE3_y] = generate_CT(eNB3_x,eNB3_y,Cell_Radius);
CT_user_list_all=[CUE1_x,CUE1_y;CUE2_x,CUE2_y;CUE3_x,CUE3_y];
[CT_user_list_all_mapped,grid_Kazan_upd] = mapping_to_grid(grid_Kazan_upd,CT_user_list_all);
grid_Kazan_for_CTs=readtable('C:\Users\Xiaomi\Documents\YASMP\ARP\codes\Points.xlsx', 'Sheet', 7);%2:-79,4:-106,5:-70
grid_Kazan_for_CTs=table2array(grid_Kazan_for_CTs);

    figure
    axis ([-650 -50 -1025 -425]);
    axis square
    grid on;
    scatter(grid_Kazan_for_CTs(:,1),grid_Kazan_for_CTs(:,2),1,'MarkerEdgeAlpha',.2,'HandleVisibility','off');
    hold on;
    scatter(CT_user_list_all_mapped(:,1),CT_user_list_all_mapped(:,2),'r','filled','d');
    hold on;
    [CT_user_list_all_mapped,~] = mapping_to_grid(grid_Kazan_for_CTs,CT_user_list_all_mapped); %remap
    scatter(CT_user_list_all_mapped(:,1),CT_user_list_all_mapped(:,2),'g','filled','d');
    title('CT SNR map');
    scatter(eNB1_x,eNB1_y,'k','filled');
    scatter(eNB2_x,eNB2_y,'k','filled');
    scatter(eNB3_x,eNB3_y,'k','filled');


% user_list=[CT_user_list_all_mapped;D2D_user_list_all_mapped_1col];
% t1='C:\Users\Xiaomi\Documents\YASMP\ARP\codes\cases\case';
% t2='\D2D_user_list_all_mapped.mat';
% save(strcat(t1,num2str(N),t2),'D2D_user_list_all_mapped');
% t2='\CT_user_list_all_mapped.mat';
% save(strcat(t1,num2str(N),t2),'CT_user_list_all_mapped'); 
% t2='\user_list.mat';
% save(strcat(t1,num2str(N),t2),'user_list');

function [user_list_mapped,grid]= mapping_to_grid (grid, user_list)
    if (isempty(user_list))
        user_list_mapped=[];
        return
    end
    if (size(user_list,2)==4) %DT Tx/Rx list w/ 4 coordinates in a row
        for i=1:size(user_list,1) 
            distances_to_grid=[];distances_to_grid_sorted=[];sort_order=[];
            for k=1:length(grid) %for Tx
                distances_to_grid(k)=dist(user_list(i,1),user_list(i,2),grid(k,1),grid(k,2));
            end     %create list with all distances from UE to each point 
            distances_to_grid=transpose(distances_to_grid);
            [distances_to_grid_sorted,sort_order]=sortrows(distances_to_grid);
            user_list_mapped(i,1)=grid(sort_order(1),1); %x
            user_list_mapped(i,2)=grid(sort_order(1),2); %y
            grid(sort_order(1),:)=[]; %current x/y point is reserved

            distances_to_grid=[];distances_to_grid_sorted=[];sort_order=[];
            for k=1:length(grid) %for Rx
                distances_to_grid(k)=dist(user_list(i,3),user_list(i,4),grid(k,1),grid(k,2));
            end %create list with all distances from UE to each point 
            distances_to_grid=transpose(distances_to_grid);
            [distances_to_grid_sorted,sort_order]=sortrows(distances_to_grid);
            user_list_mapped(i,3)=grid(sort_order(1),1); %x
            user_list_mapped(i,4)=grid(sort_order(1),2); %y
            grid(sort_order(1),:)=[];
        end
    elseif (size(user_list,2)==2) %CT list w/ 2 coordinates in a row
        for i=1:size(user_list,1) 
            distances_to_grid=[];distances_to_grid_sorted=[];sort_order=[];
            for k=1:length(grid) 
                distances_to_grid(k)=dist(user_list(i,1),user_list(i,2),grid(k,1),grid(k,2));
            end     %create list with all distances from UE to each point 
            distances_to_grid=transpose(distances_to_grid);
            [distances_to_grid_sorted,sort_order]=sortrows(distances_to_grid);
            user_list_mapped(i,1)=grid(sort_order(1),1); %x
            user_list_mapped(i,2)=grid(sort_order(1),2); %y
            grid(sort_order(1),:)=[];
        end
    end
end

function [ans] = dist (x1, y1, x2, y2) 
    ans= sqrt((abs(x1-x2))^2+(abs(y1-y2))^2);
end

function [CUE_x,CUE_y]=generate_CT(eNB_x,eNB_y,Cell_Radius)       
    UE_Dist_Min = 10; % Minimum distance of any UE (i.e. CUE or DUE) from either the BS or another UE (i.e. CUE or DUE)
    locUE = UE_Dist_Min + (Cell_Radius - UE_Dist_Min)*sqrt(rand(1,1));
    theta_= 2*pi*rand(1,1);
    CUE_x = locUE*cos(theta_) + eNB_x ;
    CUE_y = locUE*sin(theta_) + eNB_y ;
end