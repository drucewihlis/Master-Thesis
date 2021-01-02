load('C:\Users\Xiaomi\Documents\YASMP\ARP\codes\ARP_defence_codes\B1','B1');%old mPRS
load('E1','E1'); %new mAOS
load('C1','C1'); %new mPRS
figure
plot(B1(:,1),B1(:,2),'r','linewidth',2.5);%b r y m g c 
hold on
plot(E1(:,1),E1(:,2),'b','linewidth',2.5);
plot(C1(:,1),C1(:,2),'m','linewidth',2.5);
grid on
legend('old multi-cell PRS','added multi-cell AOS','upd multi-cell PRS');
xlabel('Number of D2D pairs','FontName','Arial','FontSize',14);
ylabel('Spectral Efficiency (bps/Hz)','FontName','Arial','FontSize',14);
