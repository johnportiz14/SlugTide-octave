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
noisestart1=WellInitial.noisestart1;
noiseend1=WellInitial.noiseend1;
if (nports>1)
dom_1=WellInitial.dom_1;
dom_2=WellInitial.dom_2;
dom_3=WellInitial.dom_3;
dom_4=WellInitial.dom_4;
dom_5=WellInitial.dom_5;
dom_6=WellInitial.dom_6;
dom_7=WellInitial.dom_7;
dom_8=WellInitial.dom_8;
noisestart_1=WellInitial.noisestart_1;
noiseend_1=WellInitial.noiseend_1;
noisestart_2=WellInitial.noisestart_2;
noiseend_2=WellInitial.noiseend_2;
noisestart_3=WellInitial.noisestart_3;
noiseend_3=WellInitial.noiseend_3;
noisestart_4=WellInitial.noisestart_4;
noiseend_4=WellInitial.noiseend_4;
noisestart_5=WellInitial.noisestart_5;
noiseend_5=WellInitial.noiseend_5;
noisestart_6=WellInitial.noisestart_6;
noiseend_6=WellInitial.noiseend_6;
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

        %     % if WL contains NaN (missing data from original file)   
        %     indiceNaN = isnan(Well.WL);
        %     Well.WL = Well.WL(indiceNaN~=1);
        %     Well.time = Well.time(indiceNaN~=1);

            Well.num_ports=nports;

            Well.synthfn=strvcat([TideDir tidename]);
            Well.t0_tide=[datenum(tidedate)]+tcUTC; % UTC

            Well.dt=dt; % units days

            t=Well.time;
            Well.tUTC=t+cUTC; % UTC 

            %interpolate data to preferred time interval (match data dt); change interval as necessary
            % to account for differences within original data
            xi=Well.tUTC(1):Well.dt:Well.tUTC(end);


            yi=interp1(Well.tUTC,ddata,xi);

            Well.xi=xi; % interpolated time
            Well.yi=yi; % interpolated data


            Well.date1=xi(1); % start time in UTC
            Well.date2=xi(end); % end time in UTC



            % define periods that will be eliminated in processing due to noise
            % first column is beginning of period; second is end
            Well.noise=[datenum(noisestart1) datenum(noiseend1)];

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
        %     indiceNaN = isnan(Well.WL);
        %     Well.WL = Well.WL(indiceNaN~=1);
        %     Well.time = Well.time(indiceNaN~=1);

            Well.num_ports=nports;

            Well.synthfn=strvcat([TideDir tidename]);
            Well.t0_tide=[datenum(tidedate)]+tcUTC; % UTC

            Well.dt=dt; % units days

            t=Well.time;
            Well.tUTC=t+cUTC; % UTC 

            %interpolate data to preferred time interval (match data dt); change interval as necessary
            % to account for differences within original data
            xi=Well.tUTC(1):Well.dt:Well.tUTC(end);
            yi=interp1(Well.tUTC,ddata,xi);

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


                Well.num_ports=nports;

                Well.synthfn=strvcat([TideDir tidename]);
                Well.t0_tide=[datenum(tidedate)]+tcUTC;

                Well.dt=dt; % units days

                t=Well.time;
                Well.tUTC=t+cUTC; % UTC 

                %interpolate data to preferred time interval (match data dt); change interval as necessary
                % to account for differences within original data
                xi=Well.tUTC(1):Well.dt:Well.tUTC(end);
                yi{1}=interp1(Well.tUTC,ddata_1,xi);
                yi{2}=interp1(Well.tUTC,ddata_2,xi);
                yi{3}=interp1(Well.tUTC,ddata_3,xi);
                yi{4}=interp1(Well.tUTC,ddata_4,xi);
                yi{5}=interp1(Well.tUTC,ddata_5,xi);
                yi{6}=interp1(Well.tUTC,ddata_6,xi);
                yi{7}=interp1(Well.tUTC,ddata_7,xi);
                yi{8}=interp1(Well.tUTC,ddata_8,xi);

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
        %     indiceNaN = isnan(Well.WL);
        %     Well.WL = Well.WL(indiceNaN~=1);
        %     Well.time = Well.time(indiceNaN~=1);

            Well.num_ports=nports;

            Well.synthfn=strvcat([TideDir tidename]);
            Well.t0_tide=[datenum(tidedate)]+tcUTC; % UTC

            Well.dt=dt; % units days

            t=Well.time;
            Well.tUTC=t+cUTC; % UTC 

            %interpolate data to preferred time interval (match data dt); change interval as necessary
            % to account for differences within original data
            xi=Well.tUTC(1):Well.dt:Well.tUTC(end);


            yi=interp1(Well.tUTC,ddata,xi);

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

