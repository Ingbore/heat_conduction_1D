%% Introduction
% I can't figure out dakota, so instead of that I'm going to write a matlab
% code that will find the optimal time since the step change in heat and
% the size in the step change in heat.

% Note: you should expand the size of figure 2

% Written by JWM for 3/30/2016

clear global;
figure(1)
clf
figure(2)
clf

%% Setup
%This part sets up the values that will be used in the calculations below. 

% Constants
Ts = -7.1; % Surface Temperature (C)
Qm = .04; % Heat flow (W/m^2)

k = 2; % Thermal conductivity (W/m*K)
rho = 2000; % density (kg/m^3)
c = 2000; % Specific heat capacity (J/kg*K)
kappa = k/(rho*c); % Thermal diffusivity (m^2/s)

dz = 3; % distance between z steps (m)

% Variables and Arrays
DTs = 1:0.5:6; % Array of step changes in Temperature at t = 0 (C)
jmax = length(DTs); % Number of times to iterate through the j loop

dt = dz^2/(2*kappa)/4; % maximum change in time to retain stability (s)

z = 0:dz:400; % depth array for 300 meters (m)
t = 0:dt:pi*10^7*100; % time array for 20 years (s)
imax = length(t); % number of times to iterate the time loop
T = z * Qm/k + Ts; % Initial temperature array (C)

dTdz = zeros(size(T));
dTdz(end) = .025; % Mantle heat influx (W/m^2)
Tzero = zeros(size(z));

xaxis_lower = T(1)-5; % Controls graph axes
xaxis_upper = T(end)+5; % Controls graph axes

plotcounter = 1; % Counter to plot every once in a while
tbest = 0; % variable to hold onto the best time since temperature change
DTsbest = 0; % variable to hold onto the best temperature change

load cape_thompson.dat;
capedepth = cape_thompson(:,1).'; % Depth array for Cape Thompson, transposed to be compatible with the rest of the code
capetemp = cape_thompson(:,2).'; % Temperature array for Cape Thompson, transposed to be compatible with the rest of the code

tiepoints = zeros(size(capedepth));
evaluator = 0; % Variable to determine if the current fit is the best one 
bestevaluator = 100; % Variable to track the best evaluator value (start with large value so the first evaluation overrides it)

%% Calculations
for j = 1:jmax % Goes through each possible DTs
    
    T = z * Qm/k + Ts; % Reset temperature array (C)
    T(1) = T(1) + DTs(j); % Adds step change in temperature to top (C)
    
    for i = 2:imax % Calculates each time step
    
        dTdz(1:length(T)-1) = diff(T)/dz;
        Q = -k*dTdz;
        dQdz = diff(Q)/dz;
        T(2:length(T)) = T(2:length(T)) - (1/(rho*c))*dQdz*dt;
        
        if plotcounter == 5 % Should happen once every 5 time steps
            
            plotcounter = 1; % Reset the counter
    
            figure(1)
            plot(T,z,'r','linewidth',4)
            title(['Temperature profile after ',num2str(t(i)/(pi*10^7)),' years'])
            xlabel('Temperature (C)','fontname','arial','fontsize',21)
            ylabel('Depth (m)','fontname','arial','fontsize',21)
            set(gca,'fontsize',18,'fontname','arial')
            set(gca,'YDIR','reverse')
            axis([xaxis_lower xaxis_upper z(1)-10 z(length(z))])
            pause(0.01)
        else
            plotcounter = plotcounter + 1; % Ticks up once
        end
        
        tiepoints = interp1(z,T,capedepth,'spline'); % Model temperatures at each of the Cape Thompson measured depths
        evaluator = sum(((tiepoints - capetemp).^2)); % sums up the square of the differences of the actual data and the interpolated data
        
        if evaluator < bestevaluator % Only activates when a new best has been found
            
            bestevaluator = evaluator;
            tbest = t(i);
            DTsbest = DTs(j);
            
        else
        end
        
        
    end
end
   
%% Finalize

tfinder = t(t<tbest); % creates array of proper length for following loop

T = z * Qm/k + Ts; % Reset temperature array (C)
T(1) = T(1) + DTsbest; % Adds step change in temperature to top (C)

for l = 2:length(tfinder)+1
    
    dTdz(1:length(T)-1) = diff(T)/dz;
    Q = -k*dTdz;
    dQdz = diff(Q)/dz;
    T(2:length(T)) = T(2:length(T)) - (1/(rho*c))*dQdz*dt;
    
    figure(2)
    plot(T,z,'r','linewidth',4)
    title(['Best temperature profile after ',num2str(t(l)/(pi*10^7)),' years'])
    xlabel('Temperature (C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')
    axis([xaxis_lower xaxis_upper z(1)-10 z(length(z))])
    pause(0.01)
    
end

figure (2)
hold on;
plot (capetemp,capedepth,'o')
title(['Best DT is ',num2str(DTsbest),' C and best time since change is ',num2str(tbest/(pi*10^7)),' years.'])