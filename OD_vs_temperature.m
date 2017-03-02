% Plots the optical density (absorption) for a Rb vapor cell as a
% function of temperature:
% D1 or D2 line; arbitrary laser frequency, cell length, buffer gas
% mixtures; arbitrary background dc magnetic field amplitude
%
% Optical density (OD) is defined such that Iout=Iin*exp(-OD), where Iout
% and Iin are the transmitted and incident laser intensities, respectively.
% 
% Andrew Horsley 2017
clear
units;constants;plotcolors;

%% Input parameters
purity =0.75; % abundance of 87Rb in cell (natural = 0.2783)
Dline='D2'; % choices are 'D1' or 'D2'
Bdc=0.8; % dc magnetic field in tesla
T = linspace(25,150,125)+273.15; % cell temperature
L = 2*mm; % cell length, e.g. 5*mm

bgas='N2'; %buffer gas types (as many as you want, e.g. ['Kr'; 'N2'])
Pbuffer_fill= 63*mbar; %buffer gas pressure, e.g. [10; 25]*mbar. Vector length must match bgas
Tfill=22+273.15; % temperature at which Pbuffer_fill is defined (e.g. 22+273.15 K)
n_buffer= Pbuffer_fill/kB/Tfill; % derived buffer density, in m^3

Delta = 2*pi*(266.65/2)*1e6; %laser locked to Fg=2->Fe=2/3 crossover

for tt=1:length(T)
output=Esus_func(Dline,Bdc,T(tt),L,purity,Delta,n_buffer,bgas);

% OD for each of the laser polarisations: pi (z), sigma+/- (p/m)
% OD for light polarised perpendicular to Bdc is (ODp+ODm)/2
% OD for light polarised parallel to Bdc is ODz
ODp(tt)=output.OD_p; 
ODm(tt)=output.OD_m;
ODz(tt)=output.OD_z;
end

%% Plot OD vs Temperature

fontsize=15;
linewidth=2;
legendsize=15;

xlims=[min(T-273.15),max(T-273.15)];

figure(11)

plot(T-273.15,ODp/2,'-','LineWidth',linewidth,'Color',cOrange);
hold on
plot(T-273.15,ODm/2,'-','LineWidth',linewidth,'Color',cBlue);
plot(T-273.15,ODz,'-','LineWidth',linewidth,'Color',cGreen);
hold off
set(gca,'FontSize',fontsize);
xlabel('Cell Temperature (degC)','FontSize',fontsize);
ylabel('Optical density','FontSize',fontsize);
str_ttl={sprintf('%g mm, %g mbar %s, detuning=%g MHz',L/mm,Pbuffer_fill/mbar,bgas,Delta/2/pi/MHz);...
       sprintf('Bdc=%g G',Bdc/Gauss)};
title(str_ttl,'FontWeight','Normal','FontSize',fontsize);
xlim(xlims);
leg = legend('\sigma_+','\sigma_-','\pi');
set(leg, 'Location','NorthWest','FontSize',legendsize); 

%% Plot OD spectrum
% this gives a reference to make it easier to gauge the accuracy of your
% chosen detuning

Tfixed=mean(T);
Delta_scan=linspace(-25,25,5000)*2*pi*GHz; % range of detunings to calculate OD
output=Esus_func(Dline,Bdc,Tfixed,L,purity,Delta_scan,n_buffer,bgas);

% OD for each of the laser polarisations: pi (z), sigma+/- (p/m)
% OD for light polarised perpendicular to Bdc is (ODp+ODm)/2
% OD for light polarised parallel to Bdc is ODz
ODp=output.OD_p; 
ODm=output.OD_m;
ODz=output.OD_z;

xlims=[min(Delta_scan/2/pi/GHz),max(Delta_scan/2/pi/GHz)];

figure(12)

plot(Delta_scan/2/pi/GHz,ODp,'-','LineWidth',linewidth,'Color',cOrange);
hold on
plot(Delta_scan/2/pi/GHz,ODm,'-','LineWidth',linewidth,'Color',cBlue);
plot(Delta_scan/2/pi/GHz,ODz,'-','LineWidth',linewidth,'Color',cGreen);
line([Delta Delta]/2/pi/GHz,[0 max(ODp)],'Color','k','LineWidth',linewidth);
hold off
set(gca,'FontSize',fontsize);
xlabel('Detuning from D_2 line centre (GHz)','FontSize',fontsize);
ylabel('Optical density','FontSize',fontsize);
str_ttl={sprintf('%g mm, %g degC, %g mbar %s, Bdc=%g G',L/mm, Tfixed-273.15,Pbuffer_fill/mbar,bgas,Bdc/Gauss);...
       'Black line indicates ODvsT detuning'};
title(str_ttl,'FontWeight','Normal','FontSize',fontsize);
xlim(xlims);
leg = legend('\sigma_+','\sigma_-','\pi');
set(leg, 'Location','NorthEast','FontSize',legendsize); 