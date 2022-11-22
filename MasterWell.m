clear all ;
% try
%     matlabpool
% catch
% end
% 

% Establish parameters for processing of tidal response after the well
% structure has been built
UserParam;
InitiateWells; % build the initial well input structure, inputs used to create Well structure

% %%%%%%%%%%%%%% Build Metadata Structure with waterlevel data fields
for iWell=1:NWells
    display('BuildWell...')
    iWell
    [WellArray(iWell)]=BuildWell(WellInitial(iWell));
end
Well=WellArray; % rename the full array to "Well"
save WellStruct Well;

%% %%%%%%%%%%%%Process Wells
clear;
load WellStruct;
NWells=length(Well);

for iWell=1:NWells;
    display('PlotOriginalWl...')
    iWell
    PlotOriginalWL; %plot water level
end

for iWell=1:NWells; %calculate tidal response
    display('WaterResp_general_interp...')
    iWell    
    WaterResp_general_interp;
end
% % % % % % % 
save Results Well;
% % % % 
% % % %%%Make plots
% % % 
NWells=length(Well);

for iWell=1:NWells;
    display('PlotWellResponse...')
    iWell    
    PlotWellResponse; %plot amplitude and phase
end

for iWell=1:NWells;
    display('avg_amp_phase...')
    iWell    
    avg_amp_phase; %plot port depth vs average amplitude phase
end


