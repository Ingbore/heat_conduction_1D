%%  Introduction
% Modeling of Coral Reef growth
% Written for Modeling class 2/3/16 by JWM

clear global;
figure(1)
clf

%% Initialize

x = 0:1:1000; % x array (m)
plate = -125:0.25:125; % plate depths, size of x array (m)

for j = 1:length(plate) % Creates variations in initial topography
    
    plate(j) = 7*sin(2*pi*j/100)+plate(j);
    
end

water = ones(size(x)); % water elevation size of x array (m)
coral = plate; % top of coral growth, size of x array (m)
coral(coral>water) = 0; % coral can't grow above water
growth = ones(size(length(coral))); % array of heights by which coral will grow (m)

% I used the minimum of the coefficients provided in the Galewsky paper

gmax = .01; % maximum growth rate (m/yr)
I0 = 2000; % surface light intensity (E/(s*m^2))
k = .04; % extinction coefficient (1/m)
Ik = 50; %saturating light intensity (E/(m^2*s))

dt = 1; % time step (yr)
rate = .005; % subsidence rate of plate (m/yr)
tmax = 20001; % years

xmin = x(1);
xmax = x(length(x));
ymin = plate(1)-tmax*rate;
ymax = plate(length(plate));

%% Run

for i = 1:dt:tmax
    
    water = ones(size(x));
    water = water .* 50 .* sin(2*pi*i/4000); % 50 meter oscillations, over a period of 4000 years
    
    z = water - coral; % depth of the uppermost coral (m)
    growth = gmax*tanh(I0*exp(-k*z)/Ik); % growth rate of coral over the time step (m)
    growth(coral>=water) = 0; % coral can't grow above the water
    coral = coral + growth - (dt*rate); % add growth cycle, subtract subsidence
    plate = plate - (dt*rate); % sink the plate
    
    figure(1)
    plot(x,plate,'-k','linewidth',3)
    title(['Coral reef profile on a sinking plate after ',num2str(i),' years'])
    xlabel('Distance (m)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    axis([xmin xmax ymin ymax])
    
    hold on

    plot(x,water,'--c','linewidth',1)
    plot(x,coral,'-m','linewidth',1)
    
    hold off
    pause(.00001)
    
end
    
%% Finalize