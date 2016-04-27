% this code generates a linear geotherm
% written by rsa 1/13/2016

clear all
figure(1)
clf

%% inialize
Ts = -12; % degrees C
k1 = 2.5; % W/(m-K)
k2 = 1.2;

Qm = 0.045; % W/m2

dz = 1; %meters
zmax = 800; %meters
zchange = 30; %depth of kappa change
z = 0:dz:zmax;

ksteady = k1 * ones(size(z)); %kappa array without changes

kvary = k1 * ones(size(z));
top = find(z<=30);
kvary(top) = k2; %kappa array with changes

Tsteady = Ts + (Qm ./ ksteady) .* z; %in Celsius

Tvary = Ts + (Qm ./ kvary) .* z; %gets proper slope
difference = Tvary(30) - Tsteady(30); %self explanatory
Tvary = Tvary + difference; %adds offset to entire array
Tvary(1:31) = Tvary(1:31) - difference; %subtracts offset from those that don't need it

Tzero = zeros(size(z));

%% plotting

figure(1)
plot(Tsteady,z,'r','linewidth',2)
hold on
plot(Tzero,z,'g--')
plot(Tvary,z,'b.-')


    xlabel('Temperature (°C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')

