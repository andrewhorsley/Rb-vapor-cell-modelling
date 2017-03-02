% Plots the optical density (absorption) spectrum for a Rb vapor cell:
% D1 or D2 line; arbitrary cell temperature, cell length, buffer gas
% mixtures; arbitrary background dc magnetic field amplitude
%
% Optical density (OD) is defined such that Iout=Iin*exp(-OD), where Iout
% and Iin are the transmitted and incident laser intensities, respectively.
% 
% Andrew Horsley 2017

units;constants;plotcolors;

%% Input parameters
purity =0.2783; % abundance of 87Rb in cell (natural = 0.2783)
Dline='D2'; % choices are 'D1' or 'D2'
Bdc=0.5; % dc magnetic field in tesla
T = 150+273.15; % cell temperature, e.g. 100+273.15
L = 2*mm; % cell length, e.g. 5*mm

bgas='N2'; %buffer gas types (as many as you want, e.g. ['Kr'; 'N2'])
Pbuffer_fill= 15*mbar; %buffer gas pressure, e.g. [10; 25]*mbar. Vector length must match bgas
Tfill=22+273.15; % temperature at which Pbuffer_fill is defined (e.g. 22+273.15 K)
n_buffer= Pbuffer_fill/kB/Tfill; % derived buffer density, in m^3

Delta=linspace(-25,25,5000)*2*pi*GHz; % range of detunings to calculate OD

output=Esus_func(Dline,Bdc,T,L,purity,Delta,n_buffer,bgas);

% OD for each of the laser polarisations: pi (z), sigma+/- (p/m)
% OD for light polarised perpendicular to Bdc is (ODp+ODm)/2
% OD for light polarised parallel to Bdc is ODz
ODp=output.OD_p; 
ODm=output.OD_m;
ODz=output.OD_z;

%% Plot transmission

fontsize=15;
linewidth=2;
legendsize=15;

xlims=[min(Delta/2/pi/GHz),max(Delta/2/pi/GHz)];

figure(21)

plot(Delta/2/pi/GHz,ODp,'-','LineWidth',linewidth,'Color',cOrange);
hold on
plot(Delta/2/pi/GHz,ODm,'-','LineWidth',linewidth,'Color',cBlue);
plot(Delta/2/pi/GHz,ODz,'-','LineWidth',linewidth,'Color',cGreen);
hold off
set(gca,'FontSize',fontsize);
xlabel('Detuning from D_2 line centre (GHz)','FontSize',fontsize);
ylabel('Optical density','FontSize',fontsize);
str_ttl={sprintf('%g mm, %g degC, %g mbar %s',L/mm, T-273.15,Pbuffer_fill/mbar,bgas);...
       sprintf('Bdc=%g G',Bdc/Gauss)};
title(str_ttl,'FontWeight','Normal','FontSize',fontsize);
xlim(xlims);
leg = legend('\sigma_+','\sigma_-','\pi');
set(leg, 'Location','NorthEast','FontSize',legendsize); 
