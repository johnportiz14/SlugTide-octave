function Well=LoadTides(Well)
% Load tides

UserParam;
[Ntidefiles,dum]=size(Well.synthfn);



tide=[];
ttide=[];
for i=1:Ntidefiles,
    fn=Well.synthfn(i,:);
    fn=fn(~isspace(fn));
    y=load(fn);
    y=y*1e-9;
    t=Well.t0_tide(i)+(0:length(y)-1)'*Tide_dt;
    tide=[tide; y];
    ttide=[ttide;t];
end
    
Well.tide=tide;
Well.ttide=ttide;

    
  