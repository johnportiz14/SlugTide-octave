%For wells with multiple ports (Nport > 1)
%take avg tidal response from all ports for a specific timeframe for each well

Nports=Well(iWell).num_ports;


for iport=1:Nports
      Well(iWell).T=Well(iWell).T0(:,iport);
      Well(iWell).A=Well(iWell).A0(:,iport);
      Well(iWell).Ph=Well(iWell).Ph0(:,iport);
        
        %first find phase values in timeframe for each port
        UserParam;
        [dum ind_data1]=min(abs(Well(iWell).T-datenum(PlotTime5)));
        [dum ind_data2]=min(abs(Well(iWell).T-datenum(PlotTime6)));
        ampcut=Well(iWell).A(ind_data1:ind_data2);
        phasecut=Well(iWell).Ph(ind_data1:ind_data2);

        %next find average of amplitude and phase values
        Ampc=ampcut;
        Well(iWell).Ampc=Ampc;
            % if Ampc contains NaN                                   %%
            indiceNaN = isnan(Well(iWell).Ampc);                     %%
            Well(iWell).Ampc = Well(iWell).Ampc(indiceNaN~=1);       %% 
        aAmpc=mean(Well(iWell).Ampc);
        Well(iWell).aAmpc(:,iport)=aAmpc;
        
        Phc=phasecut;
        Well(iWell).Phc=Phc;
         % if Phc contains NaN                                      %%
            indiceNaN = isnan(Well(iWell).Phc);                     %%
            Well(iWell).Phc = Well(iWell).Phc(indiceNaN~=1);        %%
        aPhc=mean(Well(iWell).Phc);
        Well(iWell).aPhc(:,iport)=aPhc;
end



%plot depth of port vs avg response
figure;clf;hold on

    plotLineMarker = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
    plotLineColor = {'b', 'k', 'r', 'g', 'm', 'c', 'k', 'b', 'r', 'g', 'k', 'r', 'g'};
    s=40;
    
for iport=1:Nports
       
        if (Nports==1)
            d=Well(iWell).dom;
        else
            d=Well(iWell).dom{iport};
        end
    

        subplot('position',[0.15,0.6,0.7,0.25])
        amp=scatter(1./Well(iWell).aAmpc(:,iport),d);
        hold on;
        title(Well(iWell).name);
        set(amp, 'MarkerEdgeColor', plotLineColor{iport}, 'Marker', plotLineMarker{iport});
        set(gca,'YDir','reverse');
        xlabel('average Amplitude (degree)');
        ylabel('Depth (m)');
       
        
        subplot('position',[0.15,0.28,0.7,0.25])
        ph=scatter(Well(iWell).aPhc(:,iport),d);
        hold on;
        set(ph, 'MarkerEdgeColor', plotLineColor{iport}, 'Marker', plotLineMarker{iport});
        set(gca,'YDir','reverse');
        xlabel('average Phase (degree)');
        ylabel('Depth (m)');
       
        
        print(gcf,'-dpdf',strcat(Well(iWell).name,'_',num2str(iport),'_avg_amp_pha.pdf'))
        print(gcf,'-dpng',strcat(Well(iWell).name,'_',num2str(iport),'_avg_amp_pha.png'))
     
    end