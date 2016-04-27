% sinusoidal forcing of a 1D thermal profile
% written by rsa for modeling class jan 2016

clear all
figure(1)
clf
figure(2)
clf
figure(3)
clf
figure(4)
clf
figure(5)
clf

%% initialize

Ts = -10; %average temperature, degrees C
DT = 15; %annual variance in temperature
kappa = 1 * 10^-6; %m2/s

dz = 0.05; %meters
zmax = 15; %meters
z=0:dz:zmax; %creates array of z
Tzero = zeros(size(z));

dtyear = 3600*24; %sec in day
Pyear = 365*3600*24; %sec in year
tyear = 0:dtyear:Pyear; %t array for year
zstaryear = sqrt(kappa*Pyear/pi);

% evaluate envelopes enclosing all temperatures

Tenv_low = Ts-DT*exp(-z/zstaryear); %lowest temperature with annual oscillation
Tenv_high = Ts+DT*exp(-z/zstaryear); %highest temperature with annual oscillation

%% run i

for i = 1:length(tyear)
    
    figure(1)
    hold on
    
    if (rem(tyear(i),dtyear*15) == 0)
        
       Tannual = Ts + DT*exp(-z./zstaryear).*sin((2*pi*tyear(i)/Pyear)-(z./zstaryear));
       plot(Tannual,z,'r','linewidth',1)
       xlabel('Temperature (°C)','fontname','arial','fontsize',21)
       ylabel('Depth (m)','fontname','arial','fontsize',21)
       set(gca,'fontsize',18,'fontname','arial')
       set(gca,'YDIR','reverse')
       axis([Ts-DT-1 Ts+DT+1 0 zmax])
       plot(Tenv_low,z,'g','linewidth',2)
       plot(Tenv_high,z,'g','linewidth',2)
       plot(Tzero,z,'k--','linewidth',2)
       
    end    
end
hold off

%% initialize ii

depths = 0:zmax/10:zmax; %array of depths throughout the range
Tzero = zeros(size(tyear)); %reset zero array to proper size

%% run ii

figure(2)
hold on
for j = 1:length(depths)
    
    Tdepth = Ts + DT*exp(-depths(j)./zstaryear).*sin((2*pi*tyear/Pyear)-(depths(j)./zstaryear));
    plot(tyear,Tdepth,'b','linewidth',1)
    xlabel('time (days)','fontname','arial','fontsize',21)
    ylabel('Temperature (°C)','fontname','arial','fontsize',21)
    plot(tyear,Tzero,'k--','linewidth',2)
    
end

%% initialize daily calculation

dT = 10; %daily temperature oscillation
dtdaily = 3600*4; %4 temperature measurements per day
tdaily = 0:dtdaily:Pyear; %t array for daily temperatures
Pdaily = 3600*24; %length of day in seconds
zstardaily = sqrt(kappa*Pdaily/pi);
Tzero = zeros(length(z)); %reset zero array to proper length

Tenv_low_daily = Tenv_low-dT*exp(-z/zstardaily); %lowest T annual+daily
Tenv_high_daily = Tenv_high+dT*exp(-z/zstardaily); %highest T annual+daily
%% run daily calculation

for l = 1:length(tdaily)
    
    figure(3)
    hold on
    
    if (rem(tdaily(l),dtdaily*65) == 0)
        
       Tannual = Ts + DT*exp(-z./zstaryear).*sin((2*pi*tdaily(l)/Pyear)-(z./zstaryear));
       Tdaily = Tannual + dT*exp(-z/zstardaily).*sin((2*pi*tdaily(l)/Pdaily)-(z./zstardaily));
       plot(Tdaily,z,'r','linewidth',1)
       xlabel('Temperature (°C)','fontname','arial','fontsize',21)
       ylabel('Depth (m)','fontname','arial','fontsize',21)
       set(gca,'fontsize',18,'fontname','arial')
       set(gca,'YDIR','reverse')
       axis([Ts-DT-1 Ts+DT+1 0 zmax])
       plot(Tenv_low_daily,z,'g','linewidth',2)
       plot(Tenv_high_daily,z,'g','linewidth',2)
       plot(Tzero,z,'k--','linewidth',2)
       
    end    
end


%% initialize movies

nframe = 0; %movie frame counter

%% run movies

for m = 1:length(tyear)
    
    figure(4)
    
    if (rem(tyear(m),dtyear*5) == 0)
        
       nframe = nframe+1;
       Tannual = Ts + DT*exp(-z./zstaryear).*sin((2*pi*tyear(m)/Pyear)-(z./zstaryear));
       plot(Tannual,z,'r','linewidth',1)
       xlabel('Temperature (°C)','fontname','arial','fontsize',21)
       ylabel('Depth (m)','fontname','arial','fontsize',21)
       set(gca,'fontsize',18,'fontname','arial')
       set(gca,'YDIR','reverse')
       axis([Ts-DT-1 Ts+DT+1 0 zmax])
       hold on
       plot(Tenv_low,z,'g','linewidth',1)
       plot(Tenv_high,z,'g','linewidth',1)
       plot(Tzero,z,'k--','linewidth',1)
       
%       M(:,nframe) = getframe(gcf);
%       pause(0.03)
       hold off
       
    end    
end

%movie2avi(M,'Tprofile_annual','fps',24)

nframe = 0; %reset movie frame counter

for n = 1:length(tdaily)
    
    figure(5)
    
    if (rem(tdaily(n),dtdaily) == 0)
       nframe = nframe+1;
       Tannual = Ts + DT*exp(-z./zstaryear).*sin((2*pi*tdaily(n)/Pyear)-(z./zstaryear));
       Tdaily = Tannual + dT*exp(-z/zstardaily).*sin((2*pi*tdaily(n)/Pdaily)-(z./zstardaily));
       plot(Tdaily,z,'r','linewidth',1)
       xlabel('Temperature (°C)','fontname','arial','fontsize',21)
       ylabel('Depth (m)','fontname','arial','fontsize',21)
       set(gca,'fontsize',18,'fontname','arial')
       set(gca,'YDIR','reverse')
       axis([Ts-DT-1 Ts+DT+1 0 zmax])
       hold on
       plot(Tenv_low_daily,z,'g','linewidth',1)
       plot(Tenv_high_daily,z,'g','linewidth',1)
       plot(Tzero,z,'k--','linewidth',1)
       
%       M(:,nframe) = getframe(gcf);
       hold off
       
    end    
end

%movie2avi(M,'Tprofile_annual+daily','fps',60)