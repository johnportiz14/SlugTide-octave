function [Well]=BuildWell(WellInitial)

UserParam;

TideDir='TideSynth/';

% call the inputs from WellIn that will be used to create the Well structure
startrow=WellInitial.startrow;
endrow=WellInitial.endrow;
filename=WellInitial.filename;
format=WellInitial.format;
Wellname=WellInitial.Wellname;
nports=WellInitial.nports;
tidename=WellInitial.tidename;
dom=WellInitial.dom;
% noisestart1=WellInitial.noisestart1;          [JPO]
% noiseend1=WellInitial.noiseend1;              [JPO]                          %%      
% [JPO] Add extra noise slots. If not detected in InitiateWell.m, sets to some arbitrary old date %%
defaultNoiseDate_start = '1900-01-01'; %default start date for noise  %%
defaultNoiseDate_end   = '1900-01-02'; %default end   date for noise  %%
if isfield(WellInitial,'noisestart1'); noisestart1=WellInitial.noisestart1;    %%
else noisestart1=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend1'); noiseend1=WellInitial.noiseend1;    %%
else noiseend1=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart2'); noisestart2=WellInitial.noisestart2;    %%
else noisestart2=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend2'); noiseend2=WellInitial.noiseend2;    %%
else noiseend2=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart3'); noisestart3=WellInitial.noisestart3;    %%
else noisestart3=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend3'); noiseend3=WellInitial.noiseend3;    %%
else noiseend3=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart4'); noisestart4=WellInitial.noisestart4;    %%
else noisestart4=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend4'); noiseend4=WellInitial.noiseend4;    %%
else noiseend4=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart5'); noisestart5=WellInitial.noisestart5;    %%
else noisestart5=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend5'); noiseend5=WellInitial.noiseend5;    %%
else noiseend5=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart6'); noisestart6=WellInitial.noisestart6;    %%
else noisestart6=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend6'); noiseend6=WellInitial.noiseend6;    %%
else noiseend6=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart7'); noisestart7=WellInitial.noisestart7;    %%
else noisestart7=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend7'); noiseend7=WellInitial.noiseend7;    %%
else noiseend7=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart8'); noisestart8=WellInitial.noisestart8;    %%
else noisestart8=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend8'); noiseend8=WellInitial.noiseend8;    %%
else noiseend8=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart9'); noisestart9=WellInitial.noisestart9;    %%
else noisestart9=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend9'); noiseend9=WellInitial.noiseend9;    %%
else noiseend9=defaultNoiseDate_end; 
end
if isfield(WellInitial,'noisestart10'); noisestart10=WellInitial.noisestart10;    %%
else noisestart10=defaultNoiseDate_start; 
end
if isfield(WellInitial,'noiseend10'); noiseend10=WellInitial.noiseend10;    %%
else noiseend10=defaultNoiseDate_end; 
end



% noisestart2=WellInitial.noisestart2;                %%
% noiseend2=WellInitial.noiseend2;                    %%

if (nports>1)
noisestart_1=WellInitial.noisestart_1;
noiseend_1=WellInitial.noiseend_1;
noisestart_2=WellInitial.noisestart_2;
noiseend_2=WellInitial.noiseend_2;
noisestart_3_1=WellInitial.noisestart_3_1;                              %%
noiseend_3_1=WellInitial.noiseend_3_1;                                  %%
noisestart_3_2=WellInitial.noisestart_3_2;                              %%
noiseend_3_2=WellInitial.noiseend_3_2;                                  %%
noisestart_3_3=WellInitial.noisestart_3_3;                              %%
noiseend_3_3=WellInitial.noiseend_3_3;                                  %%
noisestart_3_4=WellInitial.noisestart_3_4;                              %%
noiseend_3_4=WellInitial.noiseend_3_4;                                  %%
noisestart_3_5=WellInitial.noisestart_3_5;                              %%
noiseend_3_5=WellInitial.noiseend_3_5;                                  %%
noisestart_4_1=WellInitial.noisestart_4_1;                              %%
noiseend_4_1=WellInitial.noiseend_4_1;                                  %%
noisestart_4_2=WellInitial.noisestart_4_2;                              %%
noiseend_4_2=WellInitial.noiseend_4_2;                                  %%
noisestart_5=WellInitial.noisestart_5;
noiseend_5=WellInitial.noiseend_5;
noisestart_6_1=WellInitial.noisestart_6_1;                              %%
noiseend_6_1=WellInitial.noiseend_6_1;                                  %%
noisestart_6_2=WellInitial.noisestart_6_2;                              %%
noiseend_6_2=WellInitial.noiseend_6_2;                                  %%
noisestart_6_3=WellInitial.noisestart_6_3;                              %%
noiseend_6_3=WellInitial.noiseend_6_3;                                  %%
noisestart_7=WellInitial.noisestart_7;
noiseend_7=WellInitial.noiseend_7;
noisestart_8=WellInitial.noisestart_8;
noiseend_8=WellInitial.noiseend_8;
end


if format=='S' % simple format
                [Column1,Column2,Column3,Column4,Column5,Column6,Column7,Column8,Column9,Column10,Column11,Column12,Column13,Column14,Column15,Column16,Column17,Column18,Column19,Column20,Column21,Column22,Column23,Column24,Column25,Column26,Column27,Column28,Column29,Column30,Column31,Column32,Column33,Column34,Column35,Column36,Column37] = importfile(filename, startrow, endrow, format);
                atime=Column1; % time column
                adata=Column2; % Water Evelation column


            Npoints=length(adata);

            for i=1:Npoints,            % convert data from cell to double
                btime=atime(i);
                ctime=btime{1};  %curly braces make a cell array
                dtime(i)=str2double(ctime);

                bdata=adata(i);
                cdata=bdata{1};
                ddata(i)=str2double(cdata);

            end
            
            date0=datenum('1899-12-30');
            dtime=dtime+date0; % correct for difference in serial reference dates between Excel and MATLAB
            
            
            Well.name=Wellname;
            Well.time=dtime;
            Well.WL=ddata;

% %         % [JPO] 10/27/2022 (the below NaN stuff was commented out...)
            % if WL contains NaN (missing data from original file)   
            indiceNaN = isnan(Well.WL);
            Well.WL = Well.WL(indiceNaN~=1);
            Well.time = Well.time(indiceNaN~=1);
%             %%%%%

            Well.num_ports=nports;

            Well.synthfn=strvcat([TideDir tidename]);
            Well.t0_tide=[datenum(tidedate)]+tcUTC; % UTC

            Well.dt=dt; % units days

            t=Well.time;
            Well.tUTC=t+cUTC; % UTC 

            %interpolate data to preferred time interval (match data dt); change interval as necessary
            % to account for differences within original data
            xi=Well.tUTC(1):Well.dt:Well.tUTC(end);

%             yi=interp1(Well.tUTC,ddata,xi);     %ORIGINAL
            yi=interp1(Well.tUTC,Well.WL,xi);     %[JPO] Well.WL is same as ddata but with nans removed

            Well.xi=xi; % interpolated time
            Well.yi=yi; % interpolated data


            Well.date1=xi(1); % start time in UTC
            Well.date2=xi(end); % end time in UTC



            % define periods that will be eliminated in processing due to noise
            % first column is beginning of period; second is end
%             Well.noise=[datenum(noisestart1) datenum(noiseend1)];  %<-- ORIGINAL
            Well.noise=[datenum(noisestart1) datenum(noiseend1)     %%
                        datenum(noisestart2) datenum(noiseend2)     %%
                        datenum(noisestart3) datenum(noiseend3)     %%
                        datenum(noisestart4) datenum(noiseend4)     %%
                        datenum(noisestart5) datenum(noiseend5)     %%
                        datenum(noisestart6) datenum(noiseend6)     %%
                        datenum(noisestart7) datenum(noiseend7)     %%
                        datenum(noisestart8) datenum(noiseend8)     %%
                        datenum(noisestart9) datenum(noiseend9)     %%
                        datenum(noisestart10) datenum(noiseend10)]; %<-- [JPO] NEW 1/1/2022

            Well.dom=dom;



           clearvars -except iWell Well WellArray TideDir


elseif format=='A'
            [Column1,Column2,Column3,Column4,Column5,Column6,Column7,Column8,Column9,Column10,Column11,Column12,Column13,Column14,Column15,Column16,Column17,Column18,Column19,Column20,Column21,Column22,Column23,Column24,Column25,Column26,Column27,Column28,Column29,Column30,Column31,Column32,Column33,Column34,Column35,Column36,Column37] = importfile(filename, startrow, endrow, format);
            atime=Column1; % time column
            adata=Column7; % Water Evelation column


            Npoints=length(adata);

            for i=1:Npoints,            % convert data from cell to double
                btime=atime(i);
                ctime=btime{1};
                dtime(i)=str2double(ctime);

                bdata=adata(i);
                cdata=bdata{1};
                ddata(i)=str2double(cdata);

            end
            date0=datenum('1899-12-30');
            dtime=dtime+date0; % correct for difference in serial reference dates between Excel and MATLAB

            Well.name=Wellname;
            Well.time=dtime;
            Well.WL=ddata;

        %     % if WL contains NaN (missing data from original file)   
            indiceNaN = isnan(Well.WL);
            Well.WL = Well.WL(indiceNaN~=1);
            Well.time = Well.time(indiceNaN~=1);

            Well.num_ports=nports;

            Well.synthfn=strvcat([TideDir tidename]);
            Well.t0_tide=[datenum(tidedate)]+tcUTC; % UTC

            Well.dt=dt; % units days

            t=Well.time;
            Well.tUTC=t+cUTC; % UTC 

            %interpolate data to preferred time interval (match data dt); change interval as necessary
            % to account for differences within original data
            xi=Well.tUTC(1):Well.dt:Well.tUTC(end);
%             yi=interp1(Well.tUTC,ddata,xi);     %ORIGINAL
            yi=interp1(Well.tUTC,Well.WL,xi);     %[JPO] Well.WL is same as ddata but with nans removed

            Well.xi=xi; % interpolated time
            Well.yi=yi; % interpolated data


            Well.date1=xi(1); % start time in UTC
            Well.date2=xi(end); % end time in UTC



            % define periods that will be eliminated in processing due to noise
            % first column is beginning of period; second is end
            Well.noise=[datenum(noisestart1) datenum(noiseend1)];

            Well.dom=dom;



        clearvars -except iWell Well WellArray TideDir
        

        
elseif format=='B'
                [Column1,Column2,Column3,Column4,Column5,Column6,Column7,Column8,Column9,Column10,Column11,Column12,Column13,Column14,Column15,Column16,Column17,Column18,Column19,Column20,Column21,Column22,Column23,Column24,Column25,Column26,Column27,Column28,Column29,Column30,Column31,Column32,Column33,Column34,Column35,Column36,Column37] = importfile(filename, startrow, endrow, format);
                atime=Column1; % time column
                adata_1=Column5; % Water Evelation column of first port
                adata_2=Column7; % Water Evelation column of second port
                adata_3=Column9; % Water Evelation column of third port
                adata_4=Column11; % Water Evelation column of fourth port
                adata_5=Column13; % Water Evelation column of fifth port
                adata_6=Column15; % Water Evelation column of sixth port
                adata_7=Column17; % Water Evelation column of seventh port
                adata_8=Column19; % Water Evelation column of eighth port

                Npoints=length(adata_1);

                 for i=1:Npoints,            % convert data from cell to double
                    btime=atime(i);
                    ctime=btime{1};
                    dtime(i)=str2double(ctime);

                    bdata_1=adata_1(i);
                    cdata_1=bdata_1{1};
                    ddata_1(i)=str2double(cdata_1);

                    bdata_2=adata_2(i);
                    cdata_2=bdata_2{1};
                    ddata_2(i)=str2double(cdata_2);

                    bdata_3=adata_3(i);
                    cdata_3=bdata_3{1};
                    ddata_3(i)=str2double(cdata_3);

                    bdata_4=adata_4(i);
                    cdata_4=bdata_4{1};
                    ddata_4(i)=str2double(cdata_4);

                    bdata_5=adata_5(i);
                    cdata_5=bdata_5{1};
                    ddata_5(i)=str2double(cdata_5);

                    bdata_6=adata_6(i);
                    cdata_6=bdata_6{1};
                    ddata_6(i)=str2double(cdata_6);

                    bdata_7=adata_7(i);
                    cdata_7=bdata_7{1};
                    ddata_7(i)=str2double(cdata_7);

                    bdata_8=adata_8(i);
                    cdata_8=bdata_8{1};
                    ddata_8(i)=str2double(cdata_8);

                end
                date0=datenum('1899-12-30');   
                dtime=dtime+date0; % correct difference in serial reference dates between Excel and MATLAB

            %     %If repeated times in data
            %     %Remove repeated times
            %     num=find(diff(dtime)==0);
            %     dtime(num)=[];
            %     %Remove data at repeated times
            %     ddata_1(num)=[];
            %     ddata_2(num)=[];
            %     ddata_3(num)=[];
            %     ddata_4(num)=[];
            %     ddata_5(num)=[];
            %     ddata_6(num)=[];
            %     ddata_7(num)=[];
            %     ddata_8(num)=[];

                Well.name=Wellname;
                Well.time=dtime;
                Well.WL{1}=ddata_1;
                Well.WL{2}=ddata_2;
                Well.WL{3}=ddata_3;
                Well.WL{4}=ddata_4;
                Well.WL{5}=ddata_5;
                Well.WL{6}=ddata_6;
                Well.WL{7}=ddata_7;
                Well.WL{8}=ddata_8;

    % %         % [JPO] 10/27/2022 
                % if WL contains NaN (missing data from original file)   
                indiceNaN = isnan(Well.WL);
                Well.WL{1} = Well.WL{1}(indiceNaN~=1);
                Well.WL{2} = Well.WL{2}(indiceNaN~=1);
                Well.WL{3} = Well.WL{3}(indiceNaN~=1);
                Well.WL{4} = Well.WL{4}(indiceNaN~=1);
                Well.WL{5} = Well.WL{5}(indiceNaN~=1);
                Well.WL{6} = Well.WL{6}(indiceNaN~=1);
                Well.WL{7} = Well.WL{7}(indiceNaN~=1);
                Well.WL{8} = Well.WL{8}(indiceNaN~=1);
                Well.time = Well.time(indiceNaN~=1);
        %       %%%%%

                Well.num_ports=nports;

                Well.synthfn=strvcat([TideDir tidename]);
                Well.t0_tide=[datenum(tidedate)]+tcUTC;

                Well.dt=dt; % units days

                t=Well.time;
                Well.tUTC=t+cUTC; % UTC 

                %interpolate data to preferred time interval (match data dt); change interval as necessary
                % to account for differences within original data
                xi=Well.tUTC(1):Well.dt:Well.tUTC(end);
%                 yi{1}=interp1(Well.tUTC,ddata_1,xi);
%                 yi{2}=interp1(Well.tUTC,ddata_2,xi);
%                 yi{3}=interp1(Well.tUTC,ddata_3,xi);
%                 yi{4}=interp1(Well.tUTC,ddata_4,xi);
%                 yi{5}=interp1(Well.tUTC,ddata_5,xi);
%                 yi{6}=interp1(Well.tUTC,ddata_6,xi);
%                 yi{7}=interp1(Well.tUTC,ddata_7,xi);
%                 yi{8}=interp1(Well.tUTC,ddata_8,xi);
                %%[JPO] changed to allow NaN handling
                yi{1}=interp1(Well.tUTC,Well.WL{1},xi);
                yi{2}=interp1(Well.tUTC,Well.WL{2},xi);
                yi{3}=interp1(Well.tUTC,Well.WL{3},xi);
                yi{4}=interp1(Well.tUTC,Well.WL{4},xi);
                yi{5}=interp1(Well.tUTC,Well.WL{5},xi);
                yi{6}=interp1(Well.tUTC,Well.WL{6},xi);
                yi{7}=interp1(Well.tUTC,Well.WL{7},xi);
                yi{8}=interp1(Well.tUTC,Well.WL{8},xi);

                Well.xi=xi; % interpolated time
                Well.yi{1}=yi{1}; % interpolated data
                Well.yi{2}=yi{2};
                Well.yi{3}=yi{3};
                Well.yi{4}=yi{4};
                Well.yi{5}=yi{5};
                Well.yi{6}=yi{6};
                Well.yi{7}=yi{7};
                Well.yi{8}=yi{8};

                Well.date1=xi(1); % start time in UTC
                Well.date2=xi(end); % start time in UTC 

                % define periods that will be eliminated in processing due to noise
                % first column is beginning of period; second is end
                Well.noise{1}=[datenum(noisestart_1)  datenum(noiseend_1)];
                Well.noise{2}=[datenum(noisestart_2)  datenum(noiseend_2)];
                Well.noise{3}=[datenum(noisestart_3)  datenum(noiseend_3)];
                Well.noise{4}=[datenum(noisestart_4)  datenum(noiseend_4)];
                Well.noise{5}=[datenum(noisestart_5)  datenum(noiseend_5)];
                Well.noise{6}=[datenum(noisestart_6)  datenum(noiseend_6)];
                Well.noise{7}=[datenum(noisestart_7)  datenum(noiseend_7)];
                Well.noise{8}=[datenum(noisestart_8)  datenum(noiseend_8)];

                Well.dom{1}=dom_1;
                Well.dom{2}=dom_2;
                Well.dom{3}=dom_3;
                Well.dom{4}=dom_4;
                Well.dom{5}=dom_5;
                Well.dom{6}=dom_6;
                Well.dom{7}=dom_7;
                Well.dom{8}=dom_8;

  
            clearvars -except iWell Well WellArray TideDir
    
            
            
elseif format=='S' % simple format
         D=load(filename);
         dtime=D(:,1);
         ddata=D(:,2);
            
         % atime is time column; %adata is water elevation column
         Npoints=length(ddata);


            Well.name=Wellname;
            Well.time=dtime;
            Well.WL=ddata;

        %     % if WL contains NaN (missing data from original file)   
            indiceNaN = isnan(Well.WL);
            Well.WL = Well.WL(indiceNaN~=1);
            Well.time = Well.time(indiceNaN~=1);
        %   %%%%%%%%%%%%%

            Well.num_ports=nports;

            Well.synthfn=strvcat([TideDir tidename]);
            Well.t0_tide=[datenum(tidedate)]+tcUTC; % UTC

            Well.dt=dt; % units days

            t=Well.time;
            Well.tUTC=t+cUTC; % UTC 

            %interpolate data to preferred time interval (match data dt); change interval as necessary
            % to account for differences within original data
            xi=Well.tUTC(1):Well.dt:Well.tUTC(end);


%             yi=interp1(Well.tUTC,ddata,xi);     %ORIGINAL
            yi=interp1(Well.tUTC,Well.WL,xi);     %[JPO] Well.WL is same as ddata but with nans removed

            Well.xi=xi; % interpolated time
            Well.yi=yi; % interpolated data


            Well.date1=xi(1); % start time in UTC
            Well.date2=xi(end); % end time in UTC



            % define periods that will be eliminated in processing due to noise
            % first column is beginning of period; second is end
            Well.noise=[datenum(noisestart1) datenum(noiseend1)];

            Well.dom=dom;



           clearvars -except iWell Well WellArray TideDir
            
           
           
    else ['Error: unknown format']
        
end
 
                Well.tide=[];
                Well.ttide=[];
                Wellout=LoadTides(Well);
                Well=Wellout;

                %interpolate tide to match data interval
                UserParam;

                xi=Well.ttide(1):tide_dt:Well.ttide(end);
                yi=interp1(Well.ttide,Well.tide,xi);

                Well.ttide=xi;
                Well.tide=yi;
    
end

