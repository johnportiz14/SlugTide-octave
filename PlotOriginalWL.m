
UserParam;

% EQ_t=datenum(EQtime)+cUTC; %%%%UTC %%%% EQ

Nports=Well(iWell).num_ports;

for i=1:Nports,
      
    if (Nports==1)
        data=Well(iWell).WL;
    else
      data=Well(iWell).WL{i};
    end
  
   
t=Well(iWell).tUTC;
figure
plot(t,data,'-')
hold on;
% either plot the water level with the earthquake occurence, or without (select one plot command below, comment out the other) 
% plot([EQ_t EQ_t],[min(data)-0.02 max(data)+0.02],'r--') %%%% EQ
plot([min(data)-0.02 max(data)+0.02],'r--')
title(Well(iWell).name);
xlim([datenum(PlotTime1) datenum(PlotTime2)])
ylim([min(data)-0.02 max(data)+0.02])
ylabel('Water level (feet)');
datetick('KeepLimits')

print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(i),'_ori.pdf'))
print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(i),'_ori.png'))
end



