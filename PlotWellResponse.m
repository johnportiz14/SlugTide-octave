% Plot Amplitude & Phase results for a well

Nports=Well(iWell).num_ports;

 

    if Nports==1
%Plot graphs by port
        
for i=1:Nports
Well(iWell).T=Well(iWell).T0(:,i);
Well(iWell).A=Well(iWell).A0(:,i);
Well(iWell).Ph=Well(iWell).Ph0(:,i);
ff=find(Well(iWell).Ph>180);
Well(iWell).Ph(ff)=Well(iWell).Ph(ff)-360;



UserParam;
% EQ_t=datenum(EQtime)+cUTC; %%%%UTC %%%% EQ
% t_pre=EQ_t-29.5307./2; %%%% EQ
% t_post=EQ_t+29.5307./2; %%%% EQ
% kk=1:length(Well(iWell).T); % EQ
% mm=find(Well(iWell).T(kk)>t_pre & Well(iWell).T(kk)<t_post); %%%% EQ
% Well(iWell).Ph(mm)=NaN; %%%% EQ
% Well(iWell).A(mm)=NaN; %%%% EQ
% Well(iWell).s=1./(Well(iWell).A); %%%% EQ
% 
% Well(iWell).jj1=find(Well(iWell).T(kk)>t_pre-30 & Well(iWell).T(kk)<t_pre); %%%% EQ
% Well(iWell).pre_pha=nanmean(Well(iWell).Ph(Well(iWell).jj1)); %%%% EQ
% Well(iWell).pre_amp=nanmean(Well(iWell).s(Well(iWell).jj1)); %%%% EQ
% Well(iWell).jj2=find(Well(iWell).T(kk)>t_post & Well(iWell).T(kk)<t_post+30); %%%% EQ
% Well(iWell).post_pha=nanmean(Well(iWell).Ph(Well(iWell).jj2)); %%%% EQ
% Well(iWell).post_amp=nanmean(Well(iWell).s(Well(iWell).jj2)); %%%% EQ
% 
% Well(iWell).pha_change=Well(iWell).post_pha-Well(iWell).pre_pha; %%%% EQ
% Well(iWell).amp_change=Well(iWell).post_amp-Well(iWell).pre_amp; %%%% EQ
% 
% Well(iWell).jj3=find(Well(iWell).T(kk)<t_pre); %%%% EQ
% Well(iWell).background_pha=nanstd(Well(iWell).Ph(Well(iWell).jj3)); %%%% EQ
% Well(iWell).background_amp=nanstd(Well(iWell).s(Well(iWell).jj3));  %%%% EQ

% ref(iWell) = Well(iWell).date1; %%%% EQ
%Well(iWell).t=Well(iWell).T-ref(iWell); %%%% EQ
% epi(iWell)=EQ_t-ref(iWell); %%%% EQ



figure;clf;hold on
subplot('position',[0.15,0.6,0.7,0.25])

plot(Well(iWell).T,1./Well(iWell).A,'.','MarkerSize',6)
hold on;
title(Well(iWell).name);
% plot([EQ_t EQ_t],[nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2],'r--','LineWidth',0.5);hold on %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).s(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).s(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
%xlim([0 Well(iWell).t(end)+29.5307/2]);
% xlim([datenum('PlotTime3') datenum('PlotTime4')])
% ylim([nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2])
ylabel('Amplitude (1/m)');
datetick('KeepLimits')
hold off;

subplot('position',[0.15,0.28,0.7,0.25])
plot(Well(iWell).T,Well(iWell).Ph,'.','MarkerSize',6);
hold on;
% plot([EQ_t EQ_t],[-10 10]+nanmean(Well(iWell).Ph),'r--','LineWidth',0.5);hold on %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).Ph(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).Ph(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
ylabel('Phase (degree)')
hold on
% xlim([datenum('PlotTime3') datenum('PlotTime4')])
datetick('KeepLimits')
% ylim([-10 10]+nanmean(Well(iWell).Ph))

% hold off;

print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.pdf'))
print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.png'))

% end

end



    else

        for i=1:Nports
Well(iWell).T=Well(iWell).T0(:,i);
Well(iWell).A=Well(iWell).A0(:,i);
Well(iWell).Ph=Well(iWell).Ph0(:,i);
ff=find(Well(iWell).Ph>180);
Well(iWell).Ph(ff)=Well(iWell).Ph(ff)-360;

UserParam;
% EQ_t=datenum(EQtime)+cUTC; %%%%UTC  %%%% EQ
% t_pre=EQ_t-29.5307./2; %%%% EQ
% t_post=EQ_t+29.5307./2; %%%% EQ
% kk=1:length(Well(iWell).T); %%%% EQ
% mm=find(Well(iWell).T(kk)>t_pre & Well(iWell).T(kk)<t_post); %%%% EQ
% Well(iWell).Ph(mm)=NaN; %%%% EQ
% Well(iWell).A(mm)=NaN; %%%% EQ
% Well(iWell).s=1./(Well(iWell).A); %%%% EQ
% 
% Well(iWell).jj1=find(Well(iWell).T(kk)>t_pre-30 & Well(iWell).T(kk)<t_pre); %%%% EQ
% Well(iWell).pre_pha=nanmean(Well(iWell).Ph(Well(iWell).jj1)); %%%% EQ
% Well(iWell).pre_amp=nanmean(Well(iWell).s(Well(iWell).jj1)); %%%% EQ
% Well(iWell).jj2=find(Well(iWell).T(kk)>t_post & Well(iWell).T(kk)<t_post+30); %%%% EQ
% Well(iWell).post_pha=nanmean(Well(iWell).Ph(Well(iWell).jj2)); %%%% EQ
% Well(iWell).post_amp=nanmean(Well(iWell).s(Well(iWell).jj2)); %%%% EQ
% 
% Well(iWell).pha_change=Well(iWell).post_pha-Well(iWell).pre_pha; %%%% EQ
% Well(iWell).amp_change=Well(iWell).post_amp-Well(iWell).pre_amp; %%%% EQ
% 
% Well(iWell).jj3=find(Well(iWell).T(kk)<t_pre); %%%% EQ
% Well(iWell).background_pha=nanstd(Well(iWell).Ph(Well(iWell).jj3)); %%%% EQ
% Well(iWell).background_amp=nanstd(Well(iWell).s(Well(iWell).jj3)); %%%% EQ

% ref(iWell) = Well(iWell).date1; %%%% EQ
%Well(iWell).t=Well(iWell).T-ref(iWell); %%%% EQ
% epi(iWell)=EQ_t-ref(iWell); %%%% EQ



figure;clf;hold on
subplot('position',[0.15,0.6,0.7,0.25])

plot(Well(iWell).T,1./Well(iWell).A,'.','MarkerSize',6)
hold on;
title(Well(iWell).name);
% plot([EQ_t EQ_t],[nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2],'r--','LineWidth',0.5);hold on %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).s(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).s(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
%xlim([0 Well(iWell).t(end)+29.5307/2]);
% xlim([datenum('PlotTime3') datenum('PlotTime4')])
% ylim([nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2])
ylabel('Amplitude (1/m)');
datetick('KeepLimits')
hold off;

subplot('position',[0.15,0.28,0.7,0.25])
plot(Well(iWell).T,Well(iWell).Ph,'.','MarkerSize',6);
hold on;
% plot([EQ_t EQ_t],[-10 10]+nanmean(Well(iWell).Ph),'r--','LineWidth',0.5);hold on %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).Ph(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).Ph(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
ylabel('Phase (degree)')
hold on
% xlim([datenum('PlotTime3') datenum('PlotTime4')])
datetick('KeepLimits')
% ylim([-10 10]+nanmean(Well(iWell).Ph))

% hold off;

print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.pdf'))
print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.png'))
% end

        end
  

        
        
        
%Plot one graph with all ports

plotLineMarker = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
plotLineColor = {'b', 'k', 'r', 'g', 'm', 'c', 'k', 'b', 'r', 'g', 'k', 'r', 'g'};
s=20; % size of markers

figure;clf;hold on
   for iport=1:Nports

        subplot('position',[0.15,0.6,0.7,0.25])
        amp=plot(Well(iWell).T0(:,iport),1./Well(iWell).A0(:,iport));
        hold on;
        title(Well(iWell).name);
        UserParam;
        xlim([datenum(PlotTime3) datenum(PlotTime4)])
        ylabel('Amplitude (1/m)');
        datetick('KeepLimits')
        set(amp, 'Color', plotLineColor{iport}, 'Marker', plotLineMarker{iport}); 
        
        
        subplot('position',[0.15,0.28,0.7,0.25])
        ph=plot(Well(iWell).T0(:,iport),Well(iWell).Ph0(:,iport));
        hold on;
        xlim([datenum(PlotTime3) datenum(PlotTime4)])
        ylabel('Phase (degree)')
        hold on
        datetick('KeepLimits')
        set(ph, 'Color', plotLineColor{iport}, 'Marker', plotLineMarker{iport});
       

       
       print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(iport),'_amp_pha.pdf'))
       print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(iport),'_amp_pha.png'))
   end

    end
    
   