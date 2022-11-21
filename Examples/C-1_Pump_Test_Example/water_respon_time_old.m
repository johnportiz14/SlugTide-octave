function [Respon]=water_respon_time(waterdata,Volstrain,ti,dt,t_win,n_shift,per_gap)
%%%%%%%%%%%prepare for the hormonics%%%%%%%%%%%%%%
Di=waterdata;
Vol_straini=Volstrain;
%%%semidiural%%%
w_S2=2*2*pi;
w_M2=1./(12.4206/24)*2*pi;
%%diurnal%%
w_O1=1/(1.0758)*2*pi;
w_K1=1/(0.9973)*2*pi;
%%%%%%%%%%%%%%%%%%%%%%%
Sin_M2=sin(w_M2*ti)';
Cos_M2=cos(w_M2*ti)';
Sin_S2=sin(w_S2*ti)';
Cos_S2=cos(w_S2*ti)';

Sin_O1=sin(w_O1*ti)';
Cos_O1=cos(w_O1*ti)';
Sin_K1=sin(w_K1*ti)';
Cos_K1=cos(w_K1*ti)';

Sin_M2(Di==0)=0;
Cos_M2(Di==0)=0;
Sin_S2(Di==0)=0;
Cos_S2(Di==0)=0;

Sin_O1(Di==0)=0;
Cos_O1(Di==0)=0;
Sin_K1(Di==0)=0;
Cos_K1(Di==0)=0;

%%%%%%%juest test FM2. FS2%%%%%%%%%%%%%%%%%%%%%
NFFT=round(t_win/dt);

Nshift=round(NFFT*n_shift);

if length(Di)<NFFT
    Nwin=1;
else
    Nwin=fix(abs((length(Di)-NFFT))/Nshift)+1;
end

NFFT2=round(NFFT/2);

C_M=zeros(Nwin,2);
C_Mth=zeros(Nwin,2);
C_S=zeros(Nwin,2);
C_Sth=zeros(Nwin,2);
C_K=zeros(Nwin,2);
C_Kth=zeros(Nwin,2);
C_O=zeros(Nwin,2);
C_Oth=zeros(Nwin,2);
T=zeros(Nwin,1);
for i=1:Nwin
    Icut=((i-1)*Nshift+1:(i-1)*Nshift+NFFT);
    if length(Icut)>length(Di)
        Icut=Icut(1:length(Di));             
    end
    
    Dcut=Di(Icut)';
    Volcut=Vol_straini(Icut);
   
 
   
    SinMcut=Sin_M2(Icut);
    CosMcut=Cos_M2(Icut);
    SinScut=Sin_S2(Icut);
    CosScut=Cos_S2(Icut);
    
    SinKcut=Sin_K1(Icut);
    CosKcut=Cos_K1(Icut);
     SinOcut=Sin_O1(Icut);
    CosOcut=Cos_O1(Icut);
    if Nwin==1
        T(i)=0.5*(ti(1)+ti(end));
    else
        T(i)=ti((i-1)*Nshift + NFFT2);
    end

    Izero=find(abs(Dcut)<1e-5);
    if (length(Izero)>per_gap*NFFT) % exclude windows that are more than 20% gap
        C_M(i,1)=NaN;
        C_M(i,2)=NaN;
        C_Mth(i,1)=NaN;
        C_Mth(i,2)=NaN;
        
        C_S(i,1)=NaN;
        C_S(i,2)=NaN;
        C_Sth(i,1)=NaN;
        C_Sth(i,2)=NaN;
        
        C_K(i,1)=NaN;
        C_K(i,2)=NaN;
        C_Kth(i,1)=NaN;
        C_Kth(i,2)=NaN;
        
        C_O(i,1)=NaN;
        C_O(i,2)=NaN;
        C_Oth(i,1)=NaN;
        C_Oth(i,2)=NaN;
        
    else
        Dcut=(Dcut-mean(nonzeros(Dcut))).*(Dcut~=0);
        Volcut=(Volcut-mean(nonzeros(Volcut))).*(Volcut~=0);

        Dcut=detrend(Dcut).*(Dcut~=0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Dcut=nonzeros(Dcut);
%         Volcut=nonzeros(Volcut);

        Inz=(abs(Dcut)>0);
        Dcut=Dcut(Inz);
        Volcut=Volcut(Inz);

        SinMcut=nonzeros(SinMcut);
        CosMcut=nonzeros(CosMcut);
        SinScut=nonzeros(SinScut);
        CosScut=nonzeros(CosScut);
        
         SinKcut=nonzeros(SinKcut);
        CosKcut=nonzeros(CosKcut);
         SinOcut=nonzeros(SinOcut);
        CosOcut=nonzeros(CosOcut);
         
        %%%%%%%%%%%%%%

         % LSQ fit
         G=[CosMcut SinMcut CosScut SinScut CosKcut SinKcut CosOcut SinOcut ones(length(SinMcut),1)];

        C=G\Dcut;
        C_t=G\Volcut;
        C_M(i,:)=[C(1) C(2)];
        C_Mth(i,:)=[C_t(1) C_t(2)];
        C_S(i,:)=[C(3) C(4)];
        C_Sth(i,:)=[C_t(3) C_t(4)];
        C_O(i,:)=[C(7) C(8)];
        C_Oth(i,:)=[C_t(7) C_t(8)];
        C_K(i,:)=[C(5) C(6)];
        C_Kth(i,:)=[C_t(5) C_t(6)];
        
    end
end

C_M2=C_M(:,1)-sqrt(-1)*C_M(:,2);
C_M2th=C_Mth(:,1)-sqrt(-1)*C_Mth(:,2);

C_S2=C_S(:,1)-sqrt(-1)*C_S(:,2);
C_S2th=C_Sth(:,1)-sqrt(-1)*C_Sth(:,2);

C_K1=C_K(:,1)-sqrt(-1)*C_K(:,2);
C_K1th=C_Kth(:,1)-sqrt(-1)*C_Kth(:,2);

C_O1=C_O(:,1)-sqrt(-1)*C_O(:,2);
C_O1th=C_Oth(:,1)-sqrt(-1)*C_Oth(:,2);



Respon.pha_M2=angle(C_M2)-angle(C_M2th); %unit read
Respon.amp_M2=abs(C_M2)./abs(C_M2th);

Respon.pha_S2=angle(C_S2)-angle(C_S2th);
Respon.amp_S2=abs(C_S2)./abs(C_S2th);

Respon.pha_K1=angle(C_K1)-angle(C_K1th);
Respon.amp_K1=abs(C_K1)./abs(C_K1th);

Respon.pha_O1=angle(C_O1)-angle(C_O1th);
Respon.amp_O1=abs(C_O1)./abs(C_O1th);

Respon.time=T;

Respon.th_M2pha=angle(C_M2th);
Respon.water_M2pha=angle(C_M2);
end

