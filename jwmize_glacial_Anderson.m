%% Introduction
% Modeling of glacier growth
% Written for Modeling class 2/24/16 by JWM and RAS (Thanks Bob!)

clear global;
figure(1)
clf
figure(2)
clf
figure(3)
clf

%% Initialize

% spatial controls
dx = 100; % x step (m)
xmax = 20000; % Max x value (m)
dzdx = 0.07; % z step (m)
zmax = 2000; % Max z value (m)
x = 0:dx:xmax; % horizontal array (m)
z = zmax - (dzdx*x); %equation for a straight line

% time controls
dt = 0.005; % t step (yr)
tmax = 1000; % maximum t (yr)
t = 0:dt:tmax; % t array (yr)
imax = length(t);
nplots = 50; % Number of times to get a plot. Needs to divide into tmax an integer number of times
tplot = tmax/nplots; % Time step between plots

u_slide = 0.1; % basal sliding rate (m/yr)

% climate variables
dbdz = 0.01; % This is the gradient in the mass balance profile (m/yr/m)
bcap = 2;

ELA_bar = 1600; % Net accumulation is 0 at this elevation (m)
ELA_amp = 200; % Amplitude of ELA change (m)
P = 250; % Period of ELA change (yr)

% ice properties
rho_ice = 917; % density of ice (kg/m^3)
g = 9.81; % gravity (m/yr^2)
A = 2.1*10^-16; % Flow law parameter (Crazy units that cancel perfectly)

e = exp(1); % Euler's exponential

s = zeros(size(z)-1); % Slope array (unitless)
h = zeros(size(s)); % Ice thickness array, between points on x (m)
Q = zeros(size(s)); % Flow array (m/yr)


ice = zeros(size(z)); % array of ice thicknesses (m)
accumulation = zeros(size(x)); % ice growth array, based on elevation (m/yr)

j=0;
%% Run

for i = 1:imax % Loop for each time step until maximum time
    ELA = (ELA_amp * sin(2*pi*t(i)/P)) + ELA_bar;
    
    b = dbdz * (ice+z-ELA); % Mass balance profile (m)
    b = min(b,bcap); % Max accumulation rate of 2 m/yr
    
    s = diff(ice+z)/dx; % Populate the slope array

    h = ice(1:end-1) + 0.5*diff(ice);
    
    Q = (h .* u_slide)  + (A .* (rho_ice .* g .* abs(s)).^3) .* ((h .^ 5) / 5); % Populate flow array (m/yr)
    Q = [0 Q 0 ];
    
        dhdt = b - diff(Q)/dx;
        
        ice = ice + (dhdt*dt); % now update using the dt finite diff
        ice = max(ice,0);
        
        if(rem(t(i),tplot)==0)
         j = j+1;
         
        figure(1)
    plot(x/1000,z,'k','linewidth',3)
    hold on
    plot(x/1000,z+ice,'c','linewidth',3)
    
    title(['Glacial valley after ',num2str(t(i)),' years'])
    xlabel('Distance (km)','fontname','arial','fontsize',24)
    ylabel('Elevation (m)','fontname','arial','fontsize',24)
    set(gca,'fontsize',18,'fontname','arial')
    axis([0 xmax/1000 min(z) max(z) + 200])
    hold off
    pause(0.1);
    
    figure(2)
    plot(x/1000,ice,'c','linewidth',3)
    title(['Ice thickness after ',num2str(t(i)),' years'])
    xlabel('Distance (km)','fontname','arial','fontsize',24)
    ylabel('Ice thickness (m)','fontname','arial','fontsize',24)
    axis([0 xmax/1000 0 500])

    time(j) = t(i); 
    ice_vol(j) = sum(ice)*dx;

    pause(0.1)
        end

end
%% Finalize

char_vol = zeros(size(ice_vol)); % For finding characteristic time scale
char_vol(char_vol<1) = (1-(1/e))*ice_vol(length(ice_vol));

figure(3)
hold on
plot(time,ice_vol,'-c','linewidth',2)
plot(time,char_vol,'--g','linewidth',1)
title('Total Ice Volume vs Time')
xlabel('Time (yr)','fontname','arial','fontsize',21)
ylabel('Total Volume (m^3)','fontname','arial','fontsize',21)
axis([0 tmax 0 ice_vol(length(ice_vol))*1.2])
% Done!