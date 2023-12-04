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

  %(new) MC statistics
  if MC==1;
    %Standard Deviation
    Well(iWell).A_std=Well(iWell).A0_std(:,i); % std already calculated as [1/m]
    Well(iWell).Ph_std=Well(iWell).Ph0_std(:,i);
    %95% Confidence Interval (lower, upper)  [2 columns per port]
    Well(iWell).A_ci_z=Well(iWell).A0_ci_z(:,1:2,i); % ci_z already calculated as [1/m]
    Well(iWell).Ph_ci_z=Well(iWell).Ph0_ci_z(:,1:2,i);
    ff=find(Well(iWell).Ph_ci_z(:,1)>180);
    Well(iWell).Ph_ci_z(ff,1)=Well(iWell).Ph_ci_z(ff,1)-360;
    ff=find(Well(iWell).Ph_ci_z(:,2)>180);
    Well(iWell).Ph_ci_z(ff,2)=Well(iWell).Ph_ci_z(ff,2)-360;
  end
  

UserParam;

% Start plotting
figure;clf;hold on
subplot('position',[0.15,0.6,0.7,0.25])

%AMPLITUDE RESPONSE
##  plot(Well(iWell).T,1./Well(iWell).A,'.','MarkerSize',6)
if MC==0;
  plot(Well(iWell).T,1./Well(iWell).A,'.','MarkerSize',6)
else
  cap_size=0;
  % (w/ error) CI already calculated as 1/A --> do not need to reciprocal
##  errorbar(Well(iWell).T, 1./Well(iWell).A, Well(iWell).A_std, '.','MarkerSize',6)%, 'CapSize', cap_size)
  errorbar(Well(iWell).T, 1./Well(iWell).A, Well(iWell).A_std, '.')
##  errorbar(Well(iWell).T,1./Well(iWell).A, 1./Well(iWell).A-Well(iWell).A_ci_z(:,1), Well(iWell).A_ci_z(:,2)-1./Well(iWell).A, '.')
end


hold on;
title(Well(iWell).name);

% % %%%%[JPO] loop through EQs (if present)
if exist("EQtime","var")
    for j=1:size(EQtime,1)
        EQ_t=datenum(EQtime(j,:))+cUTC; %%%%UTC %%%% EQ
        t_pre=EQ_t-t_win/2; %%%% EQ
        t_post=EQ_t+t_win/2; %%%% EQ
        kk=1:length(Well(iWell).T); % EQ
        mm=find(Well(iWell).T(kk)>t_pre & Well(iWell).T(kk)<t_post); %%%% EQ
        Well(iWell).Ph(mm)=NaN; %%%% EQ
        Well(iWell).A(mm)=NaN; %%%% EQ
        %MC
        if MC == 1;
          Well(iWell).Ph_std(mm)=NaN; %%%% EQ
          Well(iWell).A_std(mm)=NaN; %%%% EQ
          Well(iWell).Ph_ci_z(mm,1:2)=NaN; %%%% EQ
          Well(iWell).A_ci_z(mm,1:2)=NaN; %%%% EQ
        endif
        Well(iWell).s=1./(Well(iWell).A); %%%% EQ
        Well(iWell).jj1=find(Well(iWell).T(kk)>t_pre-30 & Well(iWell).T(kk)<t_pre); %%%% EQ
        Well(iWell).pre_pha=nanmean(Well(iWell).Ph(Well(iWell).jj1)); %%%% EQ
        Well(iWell).pre_amp=nanmean(Well(iWell).s(Well(iWell).jj1)); %%%% EQ
        Well(iWell).jj2=find(Well(iWell).T(kk)>t_post & Well(iWell).T(kk)<t_post+30); %%%% EQ
        Well(iWell).post_pha=nanmean(Well(iWell).Ph(Well(iWell).jj2)); %%%% EQ
        Well(iWell).post_amp=nanmean(Well(iWell).s(Well(iWell).jj2)); %%%% EQ
        Well(iWell).pha_change=Well(iWell).post_pha-Well(iWell).pre_pha; %%%% EQ
        Well(iWell).amp_change=Well(iWell).post_amp-Well(iWell).pre_amp; %%%% EQ
        Well(iWell).jj3=find(Well(iWell).T(kk)<t_pre); %%%% EQ
        Well(iWell).background_pha=nanstd(Well(iWell).Ph(Well(iWell).jj3)); %%%% EQ
        Well(iWell).background_amp=nanstd(Well(iWell).s(Well(iWell).jj3));  %%%% EQ
        ref(iWell) = Well(iWell).date1; %%%% EQ
        Well(iWell).t=Well(iWell).T-ref(iWell); %%%% EQ
        epi(iWell)=EQ_t-ref(iWell); %%%% EQ
    
        plot([EQ_t EQ_t],[nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2],'r--','LineWidth',0.5);hold on %%%% EQ
        plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).s(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
        plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).s(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
    end
end
    %xlim([0 Well(iWell).t(end)+t_win/2]);
    % xlim([datenum('PlotTime3') datenum('PlotTime4')])
    % ylim([nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2])
    ylabel('Amplitude (1/m)');
    datetick('KeepLimits')
    hold off;

    %PHASE LAG
    subplot('position',[0.15,0.28,0.7,0.25])
    if MC==0;
        plot(Well(iWell).T,Well(iWell).Ph,'.','MarkerSize',6);
    else
        errorbar(Well(iWell).T, Well(iWell).Ph, Well(iWell).Ph_std, '.')
    end
    % add a zero-line if appropriate
    xl = get(gca,'XLim');
    yl = get(gca,'YLim');
    if (yl(2)>=0 && yl(1)<=0)
        hold on;
        graycolor = [0.7 0.7 0.7];
        plot(xl, [0 0], '--', 'Color', graycolor,'LineWidth',1);
        hold on;
    end
    ylabel('Phase (degree)')
    hold on;
% % %%%%[JPO] loop through EQs (if exist)
if exist("EQtime","var")
    for j=1:size(EQtime,1)
        EQ_t=datenum(EQtime(j,:))+cUTC; %%%%UTC %%%% EQ
        t_pre=EQ_t-t_win/2; %%%% EQ
        t_post=EQ_t+t_win/2; %%%% EQ
        kk=1:length(Well(iWell).T); % EQ
        mm=find(Well(iWell).T(kk)>t_pre & Well(iWell).T(kk)<t_post); %%%% EQ
        Well(iWell).Ph(mm)=NaN; %%%% EQ
        Well(iWell).A(mm)=NaN; %%%% EQ
        %MC
        if MC == 1;
          Well(iWell).Ph_std(mm)=NaN; %%%% EQ
          Well(iWell).A_std(mm)=NaN; %%%% EQ
          Well(iWell).Ph_ci_z(mm,1:2)=NaN; %%%% EQ
          Well(iWell).A_ci_z(mm,1:2)=NaN; %%%% EQ
        endif
        Well(iWell).s=1./(Well(iWell).A); %%%% EQ
        Well(iWell).jj1=find(Well(iWell).T(kk)>t_pre-30 & Well(iWell).T(kk)<t_pre); %%%% EQ
        Well(iWell).pre_pha=nanmean(Well(iWell).Ph(Well(iWell).jj1)); %%%% EQ
        Well(iWell).pre_amp=nanmean(Well(iWell).s(Well(iWell).jj1)); %%%% EQ
        Well(iWell).jj2=find(Well(iWell).T(kk)>t_post & Well(iWell).T(kk)<t_post+30); %%%% EQ
        Well(iWell).post_pha=nanmean(Well(iWell).Ph(Well(iWell).jj2)); %%%% EQ
        Well(iWell).post_amp=nanmean(Well(iWell).s(Well(iWell).jj2)); %%%% EQ
        Well(iWell).pha_change=Well(iWell).post_pha-Well(iWell).pre_pha; %%%% EQ
        Well(iWell).amp_change=Well(iWell).post_amp-Well(iWell).pre_amp; %%%% EQ
        Well(iWell).jj3=find(Well(iWell).T(kk)<t_pre); %%%% EQ
        Well(iWell).background_pha=nanstd(Well(iWell).Ph(Well(iWell).jj3)); %%%% EQ
        Well(iWell).background_amp=nanstd(Well(iWell).s(Well(iWell).jj3));  %%%% EQ
        ref(iWell) = Well(iWell).date1; %%%% EQ
        Well(iWell).t=Well(iWell).T-ref(iWell); %%%% EQ
        epi(iWell)=EQ_t-ref(iWell); %%%% EQ
        plot([EQ_t EQ_t],[-10 10]+nanmean(Well(iWell).Ph),'r--','LineWidth',0.5);hold on %%%% EQ
        plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).Ph(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
        plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).Ph(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
        hold on
    end
end
% xlim([datenum('PlotTime3') datenum('PlotTime4')])
datetick('KeepLimits')
% ylim([-10 10]+nanmean(Well(iWell).Ph))

% hold off;

print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.pdf'))
print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.png'))

%[JPO] Write amplitude response and phase shift values to .csv (time, val(or mean) )
filename = strcat(Well(iWell).name,'_',num2str(i),'_amp.csv');
twocol = [ Well(iWell).T, 1./Well(iWell).A ];
dlmwrite(filename,twocol,'precision',20);  % csvwrite() cuts off the decimal precision
filename = strcat(Well(iWell).name,'_',num2str(i),'_pha.csv');
twocol = [ Well(iWell).T, Well(iWell).Ph ];
dlmwrite(filename,twocol,'precision',20);

%(new) MC statistics (note, mean amp and pha are same as what's written above)
if MC==1;
  % Save these stats to text file (time, mean, stdev)
  filename = strcat(Well(iWell).name,'_',num2str(i),'_amp_stats.csv');
  threecol = [ Well(iWell).T, 1./Well(iWell).A, Well(iWell).A_std ];
  dlmwrite(filename,threecol,'precision',20);  % csvwrite() cuts off the decimal precision
  filename = strcat(Well(iWell).name,'_',num2str(i),'_pha_stats.csv');
  threecol = [ Well(iWell).T, Well(iWell).Ph, Well(iWell).Ph_std ];
  dlmwrite(filename,threecol,'precision',20);  % csvwrite() cuts off the decimal precision
  % Also save lower and upper 95% Confidence Interval (time, lower_CI, upper_CI)
  filename = strcat(Well(iWell).name,'_',num2str(i),'_amp_ci95.csv');
  threecol = [ Well(iWell).T, 1./Well(iWell).A , Well(iWell).A_ci_z(:,1), Well(iWell).A_ci_z(:,2)];
  dlmwrite(filename,threecol,'precision',20);  % csvwrite() cuts off the decimal precision
  filename = strcat(Well(iWell).name,'_',num2str(i),'_pha_ci95.csv');
  threecol = [ Well(iWell).T, Well(iWell).Ph, Well(iWell).Ph_ci_z(:,1), Well(iWell).Ph_ci_z(:,2)];
  dlmwrite(filename,threecol,'precision',20);  % csvwrite() cuts off the decimal precision
  
  %Print out Mean Range of Confidence intervals
  mean_A_std  = nanmean( Well(iWell).A_std);
  mean_Ph_std = nanmean( Well(iWell).Ph_std );
  mean_A_ci_z  = nanmean( [Well(iWell).A_ci_z(:,2)] - [Well(iWell).A_ci_z(:,1)] );
  mean_Ph_ci_z = nanmean( [Well(iWell).Ph_ci_z(:,2) - Well(iWell).Ph_ci_z(:,1) ]);
  printf('\n')
  printf('Time-averaged St. Dev. of Amplitude [1/m] = %e\n', mean_A_std);
  printf('Time-averaged St. Dev. of Phase     [deg] = %d\n\n', mean_Ph_std);
  printf('Time-averaged 95%% CI of Amplitude   [1/m] = %e\n', mean_A_ci_z);
  printf('Time-averaged 95%% CI of Phase       [deg] = %d\n', mean_Ph_ci_z);
end


% end

end


%% IF THERE ARE MULTIPLE PORTS  [JPO] (I have not edited this section)  !!!
    else

        for i=1:Nports
Well(iWell).T=Well(iWell).T0(:,i);
Well(iWell).A=Well(iWell).A0(:,i);
Well(iWell).Ph=Well(iWell).Ph0(:,i);
ff=find(Well(iWell).Ph>180);
Well(iWell).Ph(ff)=Well(iWell).Ph(ff)-360;

UserParam;
if exist("EQtime","var")
    for j=1:size(EQtime,1)
%         EQ_t=datenum(EQtime)+cUTC; %%%%UTC  %%%% EQ
        EQ_t=datenum(EQtime(j,:))+cUTC; %%%%UTC %%%% EQ
        t_pre=EQ_t-t_win/2; %%%% EQ
        t_post=EQ_t+t_win/2; %%%% EQ
        kk=1:length(Well(iWell).T); %%%% EQ
        mm=find(Well(iWell).T(kk)>t_pre & Well(iWell).T(kk)<t_post); %%%% EQ
        Well(iWell).Ph(mm)=NaN; %%%% EQ
        Well(iWell).A(mm)=NaN; %%%% EQ
        Well(iWell).s=1./(Well(iWell).A); %%%% EQ
        
        Well(iWell).jj1=find(Well(iWell).T(kk)>t_pre-30 & Well(iWell).T(kk)<t_pre); %%%% EQ
        Well(iWell).pre_pha=nanmean(Well(iWell).Ph(Well(iWell).jj1)); %%%% EQ
        Well(iWell).pre_amp=nanmean(Well(iWell).s(Well(iWell).jj1)); %%%% EQ
        Well(iWell).jj2=find(Well(iWell).T(kk)>t_post & Well(iWell).T(kk)<t_post+30); %%%% EQ
        Well(iWell).post_pha=nanmean(Well(iWell).Ph(Well(iWell).jj2)); %%%% EQ
        Well(iWell).post_amp=nanmean(Well(iWell).s(Well(iWell).jj2)); %%%% EQ
        
        Well(iWell).pha_change=Well(iWell).post_pha-Well(iWell).pre_pha; %%%% EQ
        Well(iWell).amp_change=Well(iWell).post_amp-Well(iWell).pre_amp; %%%% EQ
        
        Well(iWell).jj3=find(Well(iWell).T(kk)<t_pre); %%%% EQ
        Well(iWell).background_pha=nanstd(Well(iWell).Ph(Well(iWell).jj3)); %%%% EQ
        Well(iWell).background_amp=nanstd(Well(iWell).s(Well(iWell).jj3)); %%%% EQ
        
        ref(iWell) = Well(iWell).date1; %%%% EQ
        Well(iWell).t=Well(iWell).T-ref(iWell); %%%% EQ
        epi(iWell)=EQ_t-ref(iWell); %%%% EQ

        plot([EQ_t EQ_t],[-10 10]+nanmean(Well(iWell).Ph),'r--','LineWidth',0.5);hold on %%%% EQ
        plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).Ph(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
        plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).Ph(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
        hold on
    end
end


figure;clf;hold on
subplot('position',[0.15,0.6,0.7,0.25])

plot(Well(iWell).T,1./Well(iWell).A,'.','MarkerSize',6)
hold on;
title(Well(iWell).name);
% plot([EQ_t EQ_t],[nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2],'r--','LineWidth',0.5);hold on %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).s(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).s(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
% %xlim([0 Well(iWell).t(end)+t_win/2]);
% % xlim([datenum('PlotTime3') datenum('PlotTime4')])
% % ylim([nanmean(1./(Well(iWell).A))*0.8 nanmean(1./(Well(iWell).A))*1.2])
ylabel('Amplitude (1/m)');
datetick('KeepLimits')
hold off;

subplot('position',[0.15,0.28,0.7,0.25])
plot(Well(iWell).T,Well(iWell).Ph,'.','MarkerSize',6);
ylabel('Phase (degree)')
hold on;
% plot([EQ_t EQ_t],[-10 10]+nanmean(Well(iWell).Ph),'r--','LineWidth',0.5);hold on %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj1),Well(iWell).Ph(Well(iWell).jj1),'r.','MarkerSize',6) %%%% EQ
% plot(Well(iWell).T(Well(iWell).jj2),Well(iWell).Ph(Well(iWell).jj2),'r.','MarkerSize',6) %%%% EQ
hold on
% xlim([datenum('PlotTime3') datenum('PlotTime4')])
datetick('KeepLimits')
% ylim([-10 10]+nanmean(Well(iWell).Ph))

% hold off;

print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.pdf'))
print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(i),'_amp_pha.png'))
% end

        end
  

        
        
        
%% Plot one graph with all ports

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
    
   