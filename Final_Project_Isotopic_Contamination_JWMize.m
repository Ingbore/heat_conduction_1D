%% Introduction

% The goal of this code is to model the Rubidium-Strontium and
% Samarium-Neodymium isotopic systems in a body of magma.
% The system being modeled is a cooling body of basaltic magma that has
% intruded and begun to melt and incorporate granitic material above it. As
% it cools, the basalt begins to freeze out, keeping the magma at a
% constant temperature for a time and continuing to melt overlying granite
% until the molten basalt runs out. The magma begins to cool again, still 
% melting the overlying granite before finally reaching the freezing point
% of granite, at which point the systemstalls and the magma begins to
% freeze.

% Written by JWM in April 2016

clear all;
figure(1)
clf
figure(2)
clf
figure(3)
clf
figure(4)
clf
broken = 0; % flips to 1 if things get broken

%% Initialize

% Temperature
Ts = 0; % Surface Temperature (C)
T_magma = 1500; % Temperature of the intruding molten basalt (C)
T_solid_basalt = 1200; % Temperature at which basalt solidifies (C)
T_solid_granite = 700; % Temperature at which granite solidifies (C)
Qm = 0.04; % Heat flow from the mantle (W/m^2)
k_granite = 2; % Thermal conductivity of granite (W/m*K)
rho_granite = 2700; % density of granite, (kg/m^3)
c_granite = 790; % Specific heat capacity of granite (J/kg*K)
kappa_granite = k_granite/(rho_granite*c_granite); % Thermal diffusivity (m^2/s)
L_granite = 400000; % Latent heat of crystallization of granite (J/kg)
k_basalt = 1.7;
rho_basalt = 3000;
c_basalt = 840;
kappa_basalt = k_basalt/(rho_basalt*c_basalt);
L_basalt = 400000;
dz = 40; % distance between z steps (m)
zmax = 3000; % maximum depth (m)
z = 0:dz:zmax; % depth array (m)
dt = dz^2/(2*max(kappa_granite,kappa_basalt))/4; % maximum change in time to retain stability (s)
years_to_run = 20000; % number of years to run the simulation for (yr)
tmax = pi*10^7*years_to_run+dt; % length of time to run the model (s)
t = 0:dt:tmax; % time array (s)
thermal_properties = 0; % A variable that conveys what portion of the
% magma is granitic. At 0, the magma has no granite, while at 1 the magma 
% is entirely granite. This controls the specific heat, density, and other
% factors through time.
z_magma = 2500; % The location of the top of the magma
liquid_basalt = zmax-z_magma; % Volume of liquid basalt (m^3)

T = z*Qm/k_granite + Ts; % Temperature array
T(z>2500) = T_magma; % Put a 500m body of magma along the bottom of the rock column

dTdz = zeros(size(T));
dTdz(end) = 0; % Influx of heat from the bottom (W/m^2)

% Isotopic
nd_range = 50 - 5; % Most Nd concentrations range from 50 to 5 ppm
range_143 = .1232 - .1212; % Range in 143Nd natural abundance
range_144 = .2397 - .2379; % Range in 144Nd natural abuncance

granite_nd = rand * nd_range + 10; % Randomizes the concentration of Nd in the upper layer (ppm)
granite_143 = rand * range_143 + .1212; % Randomizes 143Nd in upper layer (ppm)
granite_144 = rand * range_144 + .2379; % Randomizes 144Nd in upper layer (ppm)

magma_nd = rand * nd_range + 10;
magma_143 = rand * range_143 + .1212;
magma_144 = rand * range_144 + .2379;

sr_range = 1600 - 100; % Most Sr concentrations range from 100 to 1600 ppm
range_87 = .0714 - .0694; % Range in 87Sr natural abundance (ppm)
range_86 = .0999 - .0975; % Range in 86Sr natural abundance (ppm)

granite_sr = rand * sr_range + 100; % Randomizes concentration of Sr in upper layer (ppm)
granite_87 = rand * range_87 + .0694; % Randomizes 87Sr in upper layer (ppm)
granite_86 = rand * range_86 + .0975; % Randomizes 86Sr in upper layer (ppm)

magma_sr = rand * sr_range + 100;
magma_87 = rand * range_87 + .0694;
magma_86 = rand * range_86 + .0975;

total_143_in_magma = magma_143 * rho_basalt * (zmax - z_magma) / 10^6; % Sets initial amount of 143 in magma
total_144_in_magma = magma_144 * rho_basalt *(zmax - z_magma) / 10^6;
total_87_in_magma = magma_87 * rho_basalt * (zmax - z_magma) / 10^6;
total_86_in_magma = magma_86 * rho_basalt * (zmax - z_magma) / 10^6;

%record_Nd = zeros(length(round(t(end)/dt+1)));
record_Nd(1) = total_143_in_magma/total_144_in_magma; % Start of an array that will record the 143/144Nd ratio through time
%record_Sr = zeros(length(round(t(end)/dt+1)));
record_Sr(1) = total_87_in_magma/total_86_in_magma; % To record the 87/86Sr ratio through time
%record_t = zeros(length(round(t(end)/dt+1)));
record_t(1) = 0; % Time starts at t = 0
%record_percent_solid = zeros(length(round(t(end)/dt+1)));
record_percent_liquid(1) = 100; % Magma is initially 100% liquid

% Figure parameters
xaxis_lower_1 = T(1);
xaxis_upper_1 = T_magma+100;
min_Nd = min(granite_143/granite_144,magma_143/magma_144);
max_Nd = max(granite_143/granite_144,magma_143/magma_144);
min_Sr = min(granite_87/granite_86,magma_87/magma_86);
max_Sr = max(granite_87/granite_86,magma_87/magma_86);

%% Run

for i = 1:t(end)/dt % Calculates once each time step
    
    rho_magma = thermal_properties * rho_granite + (1-thermal_properties) * rho_basalt; % overall magma density (kg/m^3)
    c_magma = thermal_properties * c_granite + (1-thermal_properties) * c_basalt; % overall magma specific heat (J/kg*K)    
    dTdz(1:length(T)-1) = diff(T)/dz;
    Q = -k_granite*dTdz; % Flux at every point (W/m^2)
    top_index = find(z>=z_magma);
    outflux = -Q(top_index(1)) * dt; % total energy outflux from the magma body each time step (J)
    if T(end)>T_solid_basalt % Everything is still hot
        z_magma = z_magma - (outflux/rho_granite/L_granite); % New top to the magma body
        thermal_properties = liquid_basalt/(zmax-z_magma); % Mixes addition of granite into the magma body
        magma_indicies = find(z>z_magma); % gets array indicies where magma is present
        magma_energy = sum(T(magma_indicies)*dz)*rho_magma*c_magma; % Calculates total thermal energy of magma (J)
        T(magma_indicies) = magma_energy/rho_magma/c_magma/(dz*length(magma_indicies)); % Magma is well mixed, so it is all the same temperature (K)
        
        total_143_in_magma = total_143_in_magma + (outflux/L_granite)*granite_143; % Addition of 143 by melting granite
        total_144_in_magma = total_144_in_magma + (outflux/L_granite)*granite_144;
        total_87_in_magma = total_87_in_magma + (outflux/L_granite)*granite_87;
        total_86_in_magma = total_86_in_magma + (outflux/L_granite)*granite_86;
        record_Nd(i+1) = total_143_in_magma/total_144_in_magma;
        record_Sr(i+1) = total_87_in_magma/total_86_in_magma;
        record_t(i+1) = i * dt / (pi*10^7); % records the time in years
        record_percent_liquid(i+1) = ((zmax-z_magma)-(500-liquid_basalt))/(zmax-z_magma) * 100;
        
        dQdz = diff(Q)/dz;
        T(2:end) = T(2:end) - (1/(rho_granite*c_granite))*dQdz*dt;
    elseif T(end)<=T_solid_basalt && liquid_basalt>0 % Once the basalt begins to freeze in the magma
        z_magma = z_magma - (outflux/rho_granite/L_granite); % New top of magma body
        thermal_properties = liquid_basalt/(zmax-z_magma); % Mixes addition of granite into magma body
        magma_indicies = find(z>z_magma); % gets array indicies where magma is present
        ideal_energy = length(magma_indicies)*1200*dz*rho_magma*c_magma; % Calculates energy necessary to have magma body at 1200 C (J)
        frozen_basalt = (ideal_energy - (sum(T(magma_indicies)*dz)*rho_magma*c_magma))/L_basalt/rho_basalt; % Volume of basalt that needs to be frozen
        if frozen_basalt<liquid_basalt % There's still plenty of basalt left
            liquid_basalt = liquid_basalt - frozen_basalt; % remove frozen basalt
            T(magma_indicies) = T_solid_basalt; % sets magma temperature at basalt freezing temperature
            
            total_143_in_magma = total_143_in_magma + (outflux/L_granite)*granite_143 - frozen_basalt*magma_143; % Addition by melting granite, subtraction by freezing basalt
            total_144_in_magma = total_144_in_magma + (outflux/L_granite)*granite_144 - frozen_basalt*magma_144;
            total_87_in_magma = total_87_in_magma + (outflux/L_granite)*granite_87 - frozen_basalt * magma_87;
            total_86_in_magma = total_86_in_magma + (outflux/L_granite)*granite_86 - frozen_basalt * magma_86;
            record_Nd(i+1) = total_143_in_magma/total_144_in_magma;
            record_Sr(i+1) = total_87_in_magma/total_86_in_magma;
            record_t(i+1) = i * dt / (pi*10^7); % records the time in years
            record_percent_liquid(i+1) = ((zmax-z_magma)-(500-liquid_basalt))/(zmax-z_magma) * 100;
            
            dQdz = diff(Q)/dz;
            T(2:end) = T(2:end) - (1/(rho_granite * c_granite))*dQdz*dt;
        else % Running out of basalt to freeze
            magma_energy = sum(T(magma_indicies)*dz)*rho_magma*c_magma + liquid_basalt*rho_basalt*L_basalt; % Energy of the magma plus the last basalt
            T(magma_indicies) = magma_energy/rho_magma/c_magma/(dz*length(magma_indicies)); % Temperature of the well mixed magma
            
            total_143_in_magma = total_143_in_magma + (outflux/L_granite)*granite_143 - liquid_basalt*magma_143; % Addition by melting granite, subtraction by freezing basalt
            total_144_in_magma = total_144_in_magma + (outflux/L_granite)*granite_144 - liquid_basalt*magma_144;
            total_87_in_magma = total_87_in_magma + (outflux/L_granite)*granite_87 - liquid_basalt * magma_87;
            total_86_in_magma = total_86_in_magma + (outflux/L_granite)*granite_86 - liquid_basalt * magma_86;
            record_Nd(i+1) = total_143_in_magma/total_144_in_magma;
            record_Sr(i+1) = total_87_in_magma/total_86_in_magma;
            record_t(i+1) = i * dt / (pi*10^7); % records the time in years
            record_percent_liquid(i+1) = ((zmax-z_magma)-(500-liquid_basalt))/(zmax-z_magma) * 100;
            
            dQdz = diff(Q)/dz;
            T(2:end) = T(2:end) - (1/(rho_granite * c_granite))*dQdz*dt;
            liquid_basalt = 0;
        end
    elseif T(end)>T_solid_granite && T(end)<=T_solid_basalt && liquid_basalt == 0 % Basalt is frozen, but granite is still liquid
        z_magma = z_magma - (outflux/rho_granite/L_granite); % New top of magma body
        thermal_properties = liquid_basalt/(zmax-z_magma); % Should be 100% granite now...
        magma_indicies = find(z>z_magma); % gets array indicies where magma is present
        magma_energy = sum(T(magma_indicies)*dz)*rho_magma*c_magma; % Calculates total thermal energy of magma (J)
        T(magma_indicies) = magma_energy/rho_magma/c_magma/(dz*length(magma_indicies)); % Temperature of the well mixed magma
        
        total_143_in_magma = total_143_in_magma + (outflux/L_granite)*granite_143; % Addition of 143 by melting granite
        total_144_in_magma = total_144_in_magma + (outflux/L_granite)*granite_144;
        total_87_in_magma = total_87_in_magma + (outflux/L_granite)*granite_87;
        total_86_in_magma = total_86_in_magma + (outflux/L_granite)*granite_86;
        record_Nd(i+1) = total_143_in_magma/total_144_in_magma;
        record_Sr(i+1) = total_87_in_magma/total_86_in_magma;
        record_t(i+1) = i * dt / (pi*10^7); % records the time in years
        record_percent_liquid(i+1) = ((zmax-z_magma)-(500-liquid_basalt))/(zmax-z_magma) * 100;
        
        dQdz = diff(Q)/dz;
        T(2:end) = T(2:end) - (1/(rho_granite * c_granite))*dQdz*dt;
    elseif T(end)<=T_solid_granite && z_magma<3000 % Granite starts to freeze
        z_magma = z_magma + outflux/L_granite/rho_granite; % Magma freezes from outside inwards
        magma_indicies = find(z>z_magma); % gets array indicies where magma is present
        T(magma_indicies) = T_solid_granite;
        
        record_Nd(i+1) = total_143_in_magma/total_144_in_magma;
        record_Sr(i+1) = total_87_in_magma/total_86_in_magma;
        record_t(i+1) = i * dt / (pi*10^7); % records the time in years
        record_percent_liquid(i+1) = ((zmax-z_magma)-(500-liquid_basalt))/(zmax-z_magma) * 100;
        
        dQdz = diff(Q)/dz;
        T(2:end) = T(2:end) - (1/(rho_granite*c_granite))*dQdz*dt;
    elseif T(end)<=T_solid_granite && z_magma>3000 % No magma left :-( just regular heat diffusion
        record_Nd(i+1) = total_143_in_magma/total_144_in_magma;
        record_Sr(i+1) = total_87_in_magma/total_86_in_magma;
        record_t(i+1) = i * dt / (pi*10^7); % records the time in years
        record_percent_liquid(i+1) = ((zmax-z_magma)-(500-liquid_basalt))/(zmax-z_magma) * 100;
        
        dQdz = diff(Q)/dz;
        T(2:end) = T(2:end) - (1/(rho_granite*c_granite))*dQdz*dt;
    else
        broken = 1;
        break
    end
    
    T_placeholder = [0 1500]; % makes a bar
    magma_height = [z_magma z_magma]; % at the height of the topmost magma

    figure(1)
    hold off
    plot(T,z,'r','linewidth',4)
    hold on
    plot(T_placeholder,magma_height,'k','linewidth',2)
    title(['Temperature profile after ',num2str(round(t(i)/(pi*10^7))),' years'])
    xlabel('Temperature (C)','fontname','arial','fontsize',12)
    ylabel('Depth (m)','fontname','arial','fontsize',12)
    set(gca,'YDIR','reverse')
    axis([xaxis_lower_1 xaxis_upper_1 z(1)-10 z(end)])
    
    pause(.01);
end 

%% Finalize
    figure(2)
    plot(record_t,record_Nd,'b','linewidth',3)
    hold on
    plot(0,record_Nd(1),'ro','linewidth',5)
    title('143/144Nd vs time')
    xlabel('Time (yr)','fontname','arial','fontsize',12)
    ylabel('143/144Nd','fontname','arial','fontsize',12)
%    axis([0 years_to_run min_Nd-.0001 max_Nd+.0001])
    
    figure(3)
    plot(record_t,record_Sr,'y','linewidth',3)
    hold on
    plot(0,record_Sr(1),'ro','linewidth',5)
    title('87/86Sr vs time')
    xlabel('Time (yr)','fontname','arial','fontsize',12)
    ylabel('87/86Sr','fontname','arial','fontsize',12)
%    axis([0 years_to_run min_Sr-.02 max_Sr+.02])
    
    figure(4)
    plot(record_t,record_percent_liquid,'k','linewidth',3)
    title('Percent of liquid in melt vs time')
    xlabel('Time (yr)','fontname','arial','fontsize',12)
    ylabel('Percent solid (%)','fontname','arial','fontsize',12)
    axis([0 years_to_run 0 100])