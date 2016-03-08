%% Introduction
% Modeling of water flow
% Written for Modeling class 3/9/16 by JWM

clear global;
figure(1)
clf
figure(2)
clf
figure(3)
clf
figure(4)
clf

%% Initialize

% spatial controls
dx = 100; % x step (m)
xmax = 20000; % Max x value (m)
dzdx = 0.07; % z step (m)
zmax = 2000; % Max z value (m)
x = 0:dx:xmax; % horizontal array (m)
z = zmax - (dzdx*x); % equation for a straight line (m)

% time controls
dt = 1; % t step (s)
tmax = 10000; % maximum t (s)
t = 0:dt:tmax; % t array (s)
imax = length(t);
nplots = 200; % Number of times to get a plot. Needs to divide into tmax an integer number of times
tplot = tmax/nplots; % Time step between plots

% climate variables
R = .02/3600; % Rainfall rate (m/s)
I =  .015/3600; % Infiltration rate (m/s)
dWdt = R-I; % This is the net rate of water accumulation (m/s)
n = 0.02; % Manning constant for a smooth sand channel (weird units that cancel)

% water properties
rho_water = 1000; % density of ice (kg/m^3)
g = 9.81; % gravity (m/s^2)

e = exp(1); % Euler's exponential

s = dzdx; % Slope
h = zeros(size(s)); % Ice thickness array, between points on x (m)
Q = zeros(size(s)); % Flow array (m/yr)

water = zeros(size(z)); % array of water depth (m)

i=0;
j=0;

%% Run

for i = 1:imax % Loop for each time step until maximum time
        
    h = water(1:end-1) + 0.5*diff(water); % Water height between points
    ubar = 1/n.*(h).^(2/3).*s.^(0.5); % Average water flow rate
    
    Q = h.*ubar; % Flow array
    Q = [0 Q Q(end)]; % Pad the flow array on both ends. The bottom end needs to drain all the water every time step
    
    dhdt = dWdt - diff(Q)/dx; % Change in water height per unit time
        
    water = water + (dhdt*dt); % now update water depth using the dt finite diff
    water = max(water,0);
        
    if(rem(t(i),tplot)==0)
        j = j+1;
        
        figure(1) % Flow depth
        plot(x/1000,water,'c','linewidth',3)
        title(['Hillslope after ',num2str(t(i)),' seconds'])
        xlabel('Distance (km)','fontname','arial','fontsize',24)
        ylabel('Depth','fontname','arial','fontsize',24)
        set(gca,'fontsize',18,'fontname','arial')
        axis([0 xmax/1000 0 0.2])
        hold off
        
    
        figure(2) % Flow discharge
        plot(x/1000,Q(2:end),'c','linewidth',3)
        title(['Water discharge after ',num2str(t(i)),' seconds'])
        xlabel('Distance (km)','fontname','arial','fontsize',24)
        ylabel('Water depth (m)','fontname','arial','fontsize',24)
        axis([0 xmax/1000 0 500])

        time(j) = t(i);
        water_vol(j) = sum(water)*dx;
        water_discharge(j) = sum(Q)*dx;

        pause(0.1)
    end

end
%% Finalize

char_vol = zeros(size(water_vol)); % For finding characteristic time scale
char_vol(char_vol<1) = (1-(1/e))*water_vol(length(water_vol));

figure(3) % Depth over time
plot(time,water_vol,'-c','linewidth',2)
title('Total Water Volume vs Time')
xlabel('Time (s)','fontname','arial','fontsize',21)
ylabel('Total Volume (m^3)','fontname','arial','fontsize',21)
axis([0 tmax 0 10000])

figure(4) % Discharge over time
plot(time,water_discharge,'-c','linewidth',2)
title('Total Water Discharge vs Time')
xlabel('Time (s)','fontname','arial','fontsize',21)
ylabel('Total Discharge (m^3/s)','fontname','arial','fontsize',21)
axis([0 tmax 0 200])

% Done!