 dbstop if error
 y5=[];
col = {'g','y','b'};
aa=cumsum(Outage_Prob_UIP(:,1));
poi1=find(Outage_Prob_UIP(:,1)~=0);
 y4=interp1(poi1(1:end-1),aa(poi1(1:end-1))/200,1:200,'Pchip');
 plot(1:200,y4,col{1},'linewidth',2.5)
 for i=2:2:20
     if mod(i/2,2)==0
        y5(i)=y4(i*10)+0.025;
     else
        y5(i)=y4(i*10)-0.025;
     end
 end
 y5(find(y5==0))=[];
 y5=[0 y5];
 hold on
 plot(0:20:200,y5,col{2},'linewidth',2.5)
 %%
 col = {'g','y','b'};
%  figure
for j=1:2
 bb=DT_SINR_Abu_cumulative_15{j};
 max_cc=max(bb);
 min_cc=min(bb);
 num=1;
 pp=[];
 for i=min_cc:0.1:max_cc
     poi=find(bb<i);
     pp(num)=length(poi)/length(bb);
     num=num+1;
 end
 plot(min_cc:0.1:max_cc,pp,col{j},'linewidth',2.5)
 hold on
end
%  hold on
%  plot(min_cc+10:0.1:max_cc+10,pp,col{2},'linewidth',2.5)