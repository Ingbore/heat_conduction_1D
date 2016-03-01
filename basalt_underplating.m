%% Introduction
% Take number 2 on heat flow
% My goal for this assignment is to create a cooling magma chamber (fig 1)
% and an incrementally assembled igneous body to see what it does to heat
% flow in an area over a period of time (initially I'm guessing 200kyr).
% 1. set up a steady state geotherm as a temperature array, then
% 2. iteratively modify it after a change in the bottom thirty kilometers.
% 3. iteratively modify it through a series of incremental emplacements of
% magma

% Written by JWM for 3/2/2016

clear global;
figure(1)
clf
figure(2)
clf
figure(3)
clf

%% Setup
%This part sets up the values that will be used in the calculations below. 

% Constants
Ts = 0; % Surface Temperature (C)
Tmagma = 1500; % Temperature of the intruding molten basalt (C)
Qm = .04; % Heat flow from the mantle (W/m^2)
k = 2; % Thermal conductivity of granite (W/m*K)
rho = 2700; % density of granite, (kg/m^3)
c = 790; % Specific heat capacity of granite (J/kg*K)
kappa = k/(rho*c); % Thermal diffusivity (m^2/s)
dz = 250; % distance between z steps (m)
zmax = 30000; % maximum depth (m)
dt = dz^2/(2*kappa)/4; % maximum change in time to retain stability (s)
tmax = pi*10^7*200000+dt; % length of time to run the model (s) - last number is number of years to run the simulation

% Variables and Arrays
z = 0:dz:zmax; % depth array (m)

t = 0:dt:tmax; % time array (s)

T_diapir = z*Qm/k+Ts; % Temperature array for magma diapir (C)
T_incremental = z*Qm/k+Ts; % Temperature array for incremental assembly (C)
T_solidus = zeros(size(T_diapir)) + 950; % Placeholder array for temperature of the granite solidus at various depths
T_diapir_max = zeros(size(T_diapir));
T_incremental_max = zeros(size(T_incremental));

T_diapir(z>20000) = Tmagma; % Put 10 km diapir of magma at 1200 C at bottom of array (C)

dTdz = zeros(size(T_diapir));
dTdz(end) = .025; % Mantle heat influx (W/m^2)

% Figure parameters
xaxis_lower_1 = T_diapir(1);
xaxis_upper_1 = Tmagma+100;
xaxis_lower_2 = T_incremental(1);
xaxis_upper_2 = Tmagma+100;

plotter = 0:100:Tmagma+100; % Placeholder array to plot against bottom of intrusion
z_magma_1 = zeros(size(plotter)) + 20000; % depth to intrusion in the diapir case (m)
z_magma_2 = zeros(size(plotter)) + 30000; % depth to intrusion in the incremental case (m)

%% Run
% This loop handles the diapir emplacement
for i = 1:t(end)/dt % Calculates once each time step
    
    dTdz(1:length(T_diapir)-1) = diff(T_diapir)/dz;
    Q = -k*dTdz;
    dQdz = diff(Q)/dz;
    T_diapir(2:end) = T_diapir(2:end) - (1/(rho*c))*dQdz*dt;
    
    hold off
    figure(1)
    plot(T_diapir,z,'r','linewidth',4)
    hold on
    plot(T_solidus,z,'y','linewidth',2)
    plot(plotter,z_magma_1,'k','linewidth',2)
    title(['Temperature profile after ',num2str(round(t(i)/(pi*10^7))),' years'])
    xlabel('Temperature (C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')
    axis([xaxis_lower_1 xaxis_upper_1 z(1)-10 z(length(z))])
    pause(0.01);
    
    T_diapir_max(T_diapir>T_diapir_max) = T_diapir(T_diapir>T_diapir_max); % Creates a maximum temperature array throughout the time interval modeled
    
end

counter = 0; % Set up a counter
% This loop handles the incremental emplacement
for j = 1:t(end)/dt % Calculates once each time step
    
    if j*dt>=counter*20000*pi*10^7 % adds basalt layer 1 km thick every 20kyr
        
        T_incremental(30000-(1000*counter) > z & z > 29000-(1000*counter)) = Tmagma; % Adds a 1 km thick basalt layer on top of the older layers once every 20000 years
        z_magma_2(z_magma_2 > 29000-(1000*counter)) = 29000-(1000*counter); % Moves the top of the intrusion to the correct depth (m)
        counter = counter + 1; % Ticks the counter up by one
        
    end
    
    dTdz(1:length(T_incremental)-1) = diff(T_incremental)/dz;
    Q = -k*dTdz;
    dQdz = diff(Q)/dz;
    T_incremental(2:end) = T_incremental(2:end) - (1/(rho*c))*dQdz*dt;
    
    hold off
    figure(2)
    plot(T_incremental,z,'r','linewidth',4)
    hold on
    plot(T_solidus,z,'y','linewidth',2)
    plot(plotter,z_magma_2,'k','linewidth',2)
    title(['Temperature profile after ',num2str(round(t(j)/(pi*10^7))),' years'])
    xlabel('Temperature (C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')
    axis([xaxis_lower_2 xaxis_upper_2 z(1) z(end)])
    pause(0.01);
    
    T_incremental_max(T_incremental>T_incremental_max) = T_incremental(T_incremental>T_incremental_max); % Creates max temperature array
    
end
%% Finalize
figure(1)
legend('Temperature gradient','Solidus of Granite','Depth to Intrusive Body')

figure(2)
legend('Temperature gradient','Solidus of Granite','Depth to Intrusive Body')

figure(3)
plot(T_diapir_max,z,'r','linewidth',2)
hold on
plot(T_incremental_max,z,'c','linewidth',2)
title('Maximum Temperature profile throughout the model time')
xlabel('Temperature (C)','fontname','arial','fontsize',21)
ylabel('Depth (m)','fontname','arial','fontsize',21)
set(gca,'fontsize',18,'fontname','arial')
set(gca,'YDIR','reverse')
axis([xaxis_lower_2 xaxis_upper_2 z(1) z(end)])
legend('Diapir','Incremental')

% Done!