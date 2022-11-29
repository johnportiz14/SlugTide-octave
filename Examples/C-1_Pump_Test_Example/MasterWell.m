clear all ;
% try
%     matlabpool
% catch
% end
% 
% %%%%%%%%%%%%%% Determine if running in Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
% If Octave, load the required packages
if isOctave
    pkg load signal
end

% Establish parameters for processing of tidal response after the well
% structure has been built
UserParam;
InitiateWells; % build the initial well input structure, inputs used to create Well structure

% %%%%%%%%%%%%%% Build Metadata Structure with waterlevel data fields
display('BuildWell...')
for iWell=1:NWells
    iWell
    [WellArray(iWell)]=BuildWell(WellInitial(iWell));
end
Well=WellArray; % rename the full array to "Well"
save WellStruct Well;

%% %%%%%%%%%%%%Process Wells
clear;
load WellStruct;
NWells=length(Well);

display('PlotOriginalWl...')
for iWell=1:NWells;
    iWell
    PlotOriginalWL; %plot water level
end

display('WaterResp_general_interp...')
for iWell=1:NWells; %calculate tidal response
    iWell    
    WaterResp_general_interp;
end
% % % % % % % 
save Results Well;
% % % % 
% % % %%%Make plots
% % % 
NWells=length(Well);

display('PlotWellResponse...')
for iWell=1:NWells;
    iWell    
    PlotWellResponse; %plot amplitude and phase
end

display('avg_amp_phase...')
for iWell=1:NWells;
    iWell    
    avg_amp_phase; %plot port depth vs average amplitude phase
end


