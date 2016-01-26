%% Introduction
% My understanding of the goal of this assignment is to 
% 1. set up a steady state geotherm as a temperature array, then
% 2. iteratively modify it after a change in surface temperature for some
% set amount of time t.

% Written by JWM for 1/27/2016

clear all;
figure(1)
clf
figure(2)
clf

%% Setup
%This part sets up the values that will be used in the calculations below. 

% Constants
Ts = -5; % Surface Temperature, C
Qm = .065; % Heat flow, W/m^2
k = 2; % Thermal conductivity, W/m*K
rho = 2700; % density, kg/m^3
c = 2000; % Specific heat capacity, J/kg*K
kappa = k/(rho*c); % Thermal diffusivity, m^2/s
dz = 1; % distance between z steps, m
DTs = 10; % Step change in Temperature at t=0, C

% Variables and Arrays
dt = dz^2/(2*kappa)/4; % maximum change in time to retain stability, s
z = 0:dz:300; % depth array for 300 meters, m
t = 0:dt:pi*10^7*20; % time array for 100 years, s
T = z*Qm/k+Ts; % Temperature array, C
T(1) = T(1) + DTs; % Surface Temperature after change, C
dTdz = zeros(size(T));
dTdz(length(T)) = .025; % Mantle heat influx, W/m^2
Tzero = zeros(size(z));
xaxis_lower = T(2)-5;
xaxis_upper = T(length(T))+5;

%% Calculations
for i = 2:length(t) % Calculates each time step
    
    dTdz(1:length(T)-1) = diff(T)/dz;
    Q = -k*dTdz;
    dQdz = diff(Q)/dz;
    T(2:length(T)) = T(2:length(T)) - (1/(rho*c))*dQdz*dt;
    pause(0.0001);
    
    figure(1)
    plot(T,z,'r','linewidth',1)
    xlabel('Temperature (C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')
    axis([xaxis_lower xaxis_upper z(1)-10 z(length(z))])
    
end

%% Part b

Tideal = z*Qm/k+Ts; % Theoretical geothermal gradient, C
Tideal = Tideal+erfc(z)*DTs; % Gradient after a step change of DTs, C

% Tideal starts as a perfect geothermal gradient
% I don't understand how the erfc gets calculated into Tideal.
% I suspect it's time-dependent, but I'm not sure how.

    figure(2)
    plot(Tideal,z,'b','linewidth',1)
    xlabel('Temperature (C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')
    axis([xaxis_lower xaxis_upper z(1)-10 z(length(z))])
    
%% Part c

% Haven't figured out part b yet.