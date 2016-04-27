%% Introduction

% This code models the growth of sand dunes from an initial flat layer of
% sand. To get a better idea of initial creation of ripples, decrease the
% time step between plots. As the code is currently written, a maximum time
% of greater than 500k impacts isn't particularly useful.

% Written by JWM for 4/6/16

clear global;
figure(1)
clf

%% Initialize

dx = .01; % Distance between bins (m)
xmax = 100; % Distance to the end (m)
x = 0:dx:xmax; % x array (m)
z = ones(size(x)) .* 10000; % Initialize z array with 10,000 sand grains in each bin, sized to match x array
z_height = z/(10000*dx); % array of heights, based on 10000 grains filling a cubic meter of sand (m)

incidence = 15; % Angle of incidence of incoming sand grains (degrees)
launch = 1; % Distance that sand grains are launched by impact (m)

counter = 1; % Initialize counter for when to plot

%% Run

for i = 1:500000 % Pseudo-time loop calculates once for each sand impact

    min_intercept = z_height(1); % Minimum intercept value for impacting grains (m)
    z_intercept = x*tan(incidence*pi/180) + z_height; % Calculates intercept needed for impacting grain to hit each point (m)
    max_intercept = max(z_intercept); % Picks highest intercept needed for impacting grains to have a change to hit every point that can be hit (m)
    intercept_range = max_intercept - min_intercept; % Finds the range in possible intercept values (m)
    
    intercept = rand * intercept_range + min_intercept; % Picks a random intercept from the range of possible intercepts (m)
    impactor = intercept - x * tan(incidence*pi/180); % Calculates trajectory of impacting sand grain
    potential_impacts = find(z_height>impactor); % Lists all array indicies where the impacting grain can hit the sand
    impact_site = potential_impacts(1); % First potential impact site is the actual impact site
    landing_site = impact_site + launch/dx; % Array index of site where sand from impact lands
    
    if impact_site == 1 % Set of if statements to wrap the removal of sand
        front = length(z);
        back = impact_site + 1;
    elseif impact_site == length(z) % If impact happens at the end of the array
        front = impact_site - 1;
        back = 1;
    else % 
        front = impact_site - 1;
        back = impact_site + 1;
    end
    
    z(impact_site) = z(impact_site) - 12; % Removes 12 sand grains from the impact site
    z(front) = z(front) - 4; % Removes 4 sand grains from the spot in front of the impact
    z(back) = z(back) -4; % Removes 4 sand grains from the spot behind the impact
    
    if landing_site == length(z) % Set of if statements to wrap the addition of sand
        front = landing_site - 1;
        back = 1;
    elseif landing_site == length(z) + 1
        front = length(z);
        landing_site = 1;
        back = landing_site + 1;
    elseif landing_site > length(z) + 1
        landing_site = landing_site - 10001;
        front = landing_site - 1;
        back = landing_site + 1;
    else
        front = landing_site - 1;
        back = landing_site + 1;
    end
    
    z(landing_site) = z(landing_site) + 10; % Adds 10 grains of sand to the landing site
    z(front) = z(front) + 5; % Adds 5 grains of sand to the spot in front of the landing site
    z(back) = z(back) + 5; % Adds 5 grains of sand to the spot behind the landing site
    
    z_height = z/(10000*dx); % Recalculate z_height array for next loop iteration
    
    if counter == 500 % Loop to plot once every 500 time steps
        
        plot(x,z_height,'k','linewidth',2)
        title(['Temperature profile after ',num2str(i),' impacts'])
        xlabel('Distance (m)','fontname','arial','fontsize',21)
        ylabel('Height (m)','fontname','arial','fontsize',21)
        set(gca,'fontsize',18,'fontname','arial')
        axis([0 xmax 50 150])
        
        pause(.001);
        counter = 1; % Reset plotting counter
    else
        counter = counter + 1; % Ticks up counter once
    end
    
end

%% Finalize