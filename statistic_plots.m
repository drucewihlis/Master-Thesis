   clc;close all;clear;format compact;
   N=30;
   
   SINR_C_OFDMPl_PC_all=[];
   SINR_D_OFDMPl_PC_all=[];
   SINR_C_AOS_all=[];
   SINR_D_AOS_all=[];
   for n=1:N
       t1='C:\Users\Xiaomi\Documents\YASMP\ARP\codes\cases\case';
       t2=num2str(n);
       t3='\SINR_C_AOS_dB.mat';
       t=strcat(t1,t2,t3); 
       load(t,'SINR_C_AOS_dB'); 
       x=strcat('SINR_C_AOS_dB',num2str(n));
       eval([x '=SINR_C_AOS_dB;']);
       
       t3='\SINR_D_AOS_dB.mat';
       t=strcat(t1,t2,t3); 
       load(t,'SINR_D_AOS_dB'); 
       x=strcat('SINR_D_AOS_dB',num2str(n));
       eval([x '=SINR_D_AOS_dB;']);
       
        mult=1;
        for k=1:length(SINR_D_AOS_dB)
            mult=mult * ( 1 + db2pow(SINR_D_AOS_dB(k)) );
        end
        for k=1:length(SINR_C_AOS_dB)
            mult=mult * ( 1 + db2pow(SINR_C_AOS_dB(k)) );
        end
        SE_AOS(n)=log2(mult);
        qtity_of_pairs_AOS(n)=size(SINR_D_AOS_dB,1);
    
       t3='\SINR_C_OFDMPl_PC_dB.mat';
       t=strcat(t1,t2,t3); 
       load(t,'SINR_C_OFDMPl_PC_dB'); 
       x=strcat('SINR_C_OFDMPl_PC_dB',num2str(n));
       eval([x '=SINR_C_OFDMPl_PC_dB;']);
       
       t3='\SINR_D_OFDMPl_PC_dB.mat';
       t=strcat(t1,t2,t3); 
       load(t,'SINR_D_OFDMPl_PC_dB'); 
       x=strcat('SINR_D_OFDMPl_PC_dB',num2str(n));
       eval([x '=SINR_D_OFDMPl_PC_dB;']);
       
       mult=1;
       for k=1:length(SINR_D_OFDMPl_PC_dB)
           mult=mult * ( 1 + db2pow(SINR_D_OFDMPl_PC_dB(k)) );
       end
       for k=1:length(SINR_C_OFDMPl_PC_dB)
           mult=mult * ( 1 + db2pow(SINR_C_OFDMPl_PC_dB(k)) );
       end
       SE_OFDMPl_PC(n)=log2(mult);
       qtity_of_pairs_OFDMPl_PC(n)=size(SINR_D_OFDMPl_PC_dB,1);
        
       SINR_C_OFDMPl_PC_all=[SINR_C_OFDMPl_PC_all;SINR_C_OFDMPl_PC_dB];
       SINR_D_OFDMPl_PC_all=[SINR_D_OFDMPl_PC_all;SINR_D_OFDMPl_PC_dB];
       
       SINR_C_AOS_all=[SINR_C_AOS_all;SINR_C_AOS_dB];
       SINR_D_AOS_all=[SINR_D_AOS_all;SINR_D_AOS_dB];
   end

SE_OFDMPl_PC=transpose(SE_OFDMPl_PC);
SE_AOS=transpose(SE_AOS);
qtity_of_pairs_OFDMPl_PC=transpose(qtity_of_pairs_OFDMPl_PC);
qtity_of_pairs_AOS=transpose(qtity_of_pairs_AOS);

%CDF DT
figure
cdfplot(SINR_D_OFDMPl_PC_all);
xlim([0 50]);
hold on

cdfplot(SINR_D_AOS_all);
grid on
legend('OFDMPlanning PC','multi-cell AOS');
xlabel('DT SINR (dB)','FontName','Arial','FontSize',14);
ylabel('CDF','FontName','Arial','FontSize',14);

%CDF CT
figure
cdfplot(SINR_C_OFDMPl_PC_all);
xlim([0 50]);
hold on

cdfplot(SINR_C_AOS_all);
grid on
legend('OFDMPlanning PC','multi-cell AOS');
xlabel('CT SINR (dB)','FontName','Arial','FontSize',14);
ylabel('CDF','FontName','Arial','FontSize',14);

%SE
A=[qtity_of_pairs_AOS,SE_AOS]; 
B=[qtity_of_pairs_OFDMPl_PC,SE_OFDMPl_PC];
A=sortrows(A); %from smallest # of pairs to highest
[Alink,aa] = findgroups(A(:,1)); %aa - different values of pairs, Alink - links from aa to A
A1 = [ aa, splitapply(@mean,A(:,2),Alink)]; %find mean value for each # of pairs
B=sortrows(B);
[Blink,bb] = findgroups(B(:,1));
B1 = [ bb, splitapply(@mean,B(:,2),Blink)];

figure
plot(A1(:,1),A1(:,2),'y-o','linewidth',2.5); %b r y m g c 
hold on
plot(B1(:,1),B1(:,2),'g-o','linewidth',2.5);
grid on
legend('multi-cell AOS','OFDMPlanning PC');
xlabel('Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);


    
