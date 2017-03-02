% Produces 2D plots of microwave magnetic field sensitivity in a vapor cell
% as a function of buffer gas filling and vapor cell temperature. Both the
% photon-shot-noise and atomic-projection-noise limited sensitivities are
% given. The script also provides overlays of the corresponding spatial
% resolution, and plots of the optical density and T1 time
%
% Andrew Horsley, 2017

clear
constants;
units;

%% Input parameters

%%%%% general details %%%%%
InputParam.rb87_abundance=0.98; % abundance of 87Rb in cell (natural = 0.2783)
InputParam.Dline='D2'; % choices are 'D1' or 'D2'
InputParam.Bdc=0.8; % dc magnetic field in tesla
InputParam.DutyCycle = 500; % long-term duty cycle, shots/second

%%%%% laser properties %%%%%
InputParam.I0=30*mW/cm^2; % incident laser intensity
InputParam.Tprobe = 0.3*mus; %s, probe pulse length
% detuning of laser from centre of D1/D2 line. Note that this should be set
% to the desired spectral feature BEFORE the buffer gas shift. In order to
% account for the change in buffer gas shift as cell filling parameters are
% varied, the buffer gas shift is automatically compensated for in 
% sensitivity_Bmw_func. Best to check using OD_spectra with Pfill=0
% -2*pi*2.395*GHz for Bdc=0, addressing F=2 state
% 2*pi*6.646*GHz for Bdc=0.8T, addressing optical T8-pi transition
InputParam.detuning = 2*pi*6.646*GHz; 

%%%%% population difference %%%%%
% between states coupled by mw, as fraction of total population in all states. 
% For Bdc=0 and 100% optically pumped population equally populating each 
% of the F=1 states, TransPopDiff=1/3. 
% For Bdc=1T and complete depopulation of one of the coupled states,
% with equal redistribution to all 7 other states, TransPopDiff = 9/64. 
InputParam.TransPopDiff=9/64; 

%%%%% camera properties %%%%%
InputParam.QE = 0.35; % quantum efficiency of camera
InputParam.gain=1; % camera gain. default set to 1

%%%%% temperature scan %%%%%
T=linspace(80,170,50)+273.15; % cell temperature, in Kelvin
InputParam.Tfill = 22+ 273.15; %standard temperature (eg buffer gas pressure is 10mbar @ T0) - 0degC from BOlsen thesis

%%%%% Buffer pressure scan %%%%%
InputParam.bgas=[sprintf('Kr'); sprintf('N2')]; % buffer gas types (as many as you want)
Ptotal = linspace(0,300,50)*mbar;
P1_frac = 0.75; P2_frac = 1-P1_frac;
P0_1 = P1_frac*Ptotal; P0_2 = P2_frac*Ptotal;
Pfill=[P0_1; P0_2];

%%%%% Cell dimensions %%%%%
InputParam.cell_shape = 'rectangle'; % 'sphere', 'cylinder', rectangle'
if strcmp(InputParam.cell_shape, 'sphere')
    cell_radius = 0.5*mm; % 1/2 * optical path length
    InputParam.cell_dimensions = cell_radius;
elseif strcmp(InputParam.cell_shape, 'cylinder')
    cell_radius = 0.5*mm; 
    cell_thickness = 37.5*mm; % optical path length
    InputParam.cell_dimensions = [cell_radius, cell_thickness];
elseif strcmp(InputParam.cell_shape, 'rectangle')
    a = 6*mm; % transverse dimensions
    b = 6*mm; % transverse dimensions
    c = 200*mum; % cell thickness (optical path length)
    InputParam.cell_dimensions = [a, b, c];
end

%%   
%prepare arrays
T1_array=zeros(size(T,1),size(Pfill,1));
OD87_array=zeros(size(T,1),size(Pfill,1));
delta_array=zeros(size(T,1),size(Pfill,1));
Sensitivity=zeros(size(T,1),size(Pfill,1));
dBphoton=zeros(size(T,1),size(Pfill,1));
dBatom=zeros(size(T,1),size(Pfill,1));
spatial_res=zeros(size(T,1),size(Pfill,1));
Omega_min=zeros(size(T,1),size(Pfill,1));
Omega_min_opt=zeros(size(T,1),size(Pfill,1));
Nat=zeros(size(T,1),size(Pfill,1));
Ncounts=zeros(size(T,1),size(Pfill,1));
delta_Ncounts=zeros(size(T,1),size(Pfill,1));
OD=zeros(size(T,1),size(Pfill,1));
OD87=zeros(size(T,1),size(Pfill,1));

for k1=1:size(T,2)
for k2=1:size(Ptotal,2)

    temp=sensitivity_Bmw_func(T(k1),Pfill(:,k2),InputParam);
    
    Sensitivity(k2,k1)=temp.dBphoton;
    dBphoton(k2,k1)=temp.dBphoton_opt;
    dBatom(k2,k1)=temp.dBatom;
    
    spatial_res(k2,k1) =temp.diffusion_length;
    Omega_min(k2,k1)=temp.Omega_min;
    Omega_min_opt(k2,k1)=temp.Omega_min_opt;
    Nat(k2,k1)=temp.Nat;
    Ncounts(k2,k1)=temp.Ncounts;
    delta_Ncounts(k2,k1)=temp.delta_Ncounts;
    OD(k2,k1)=temp.OD;
    OD87(k2,k1)=temp.OD87;
    T1_array(k2,k1) = temp.T1;
end
end


%% Plotting: photon shot noise and atomic projection noise
% contour lines show T1-limited spatial resolution (diffusion distance
% during T1)
fontsize=14;

figure(13)
colormap(jet);

h_vector=T-273.15;
v_vector=Ptotal/mbar;

subplot(1,2,1)
imagesc(h_vector,v_vector,dBatom/pT); axis tight
c_bar=colorbar;
set(c_bar,'FontSize',fontsize);
caxis([0, 5])

hold on
contourlines=linspace(0, 100,21);
[C h]=contour(h_vector,v_vector,spatial_res/mum,contourlines,'LineColor','w');
hold off
clabel(C,h,'Color','w')
set(gca,'YDir','normal','XDir','normal')
set(gca,'FontSize',fontsize);

xlabel('Temperature (^{\circ}C)');ylabel('Total P_{fill} (mbar)');
title({'a) \deltaB_{atom}'...
    '(contour lines show T1-limited spatial resolution)'},...
    'FontSize',fontsize,'FontWeight','Normal');
annotation('textbox',...
    [0.405 0.05 0.2 0.065],...
    'String',{'$\delta B_{atom}$';'(pT Hz$^{-1/2}$)'},...
    'FontSize',fontsize-1,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'interpreter', 'latex',...
    'LineStyle','none',...
    'EdgeColor',[1 1 1]);

subplot(1,2,2)
imagesc(h_vector,v_vector,dBphoton/nT); axis tight
c_bar=colorbar;
set(c_bar,'FontSize',fontsize);
caxis([0, 100])

hold on
contourlines=linspace(0, 100,21);
[C h]=contour(h_vector,v_vector,spatial_res/mum,contourlines,'LineColor','w');
hold off
clabel(C,h,'Color','w')
set(gca,'YDir','normal','XDir','normal')
set(gca,'FontSize',fontsize);

xlabel('Temperature (^{\circ}C)');ylabel('Total P_{fill} (mbar)');
title({'b) \deltaB_{photon}'...
    '(contour lines show T1-limited spatial resolution)'},...
    'FontSize',fontsize,'FontWeight','Normal');

annotation('textbox',...
    [0.84 0.05 0.2 0.065],...
    'String',{'$\delta B_{photon}$';'(nT Hz$^{-1/2}$)'},...
    'FontSize',fontsize-1,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'interpreter', 'latex',...
    'LineStyle','none',...
    'EdgeColor',[1 1 1]);


%%
set(gcf,'PaperPosition', [-2 0 32 12]);
set(gcf,'PaperSize', [28 12]);

% print -dpdf buffer_Tstem_vs_Batom_Bphoton_87
% open buffer_Tstem_vs_Batom_Bphoton_87.pdf

%% Plotting: Optical density (OD87) and T1 time

figure(15)
colormap(jet);
imagesc(h_vector,v_vector,OD87); axis tight
c_bar=colorbar;
set(c_bar,'FontSize',fontsize);
caxis([0, 5])

hold on
contourlines=linspace(0, 50,26);
[C h]=contour(h_vector,v_vector,T1_array/mus,contourlines,'LineColor','w');
hold off
clabel(C,h,'Color','w')
set(gca,'YDir','normal','XDir','normal')
set(gca,'FontSize',fontsize);

xlabel('Temperature (^{\circ}C)');Hylabel=ylabel('Total P_{fill} (mbar)');
title('OD_{87}, with T1 contour lines','FontSize',fontsize,'FontWeight','Normal');

annotation('textbox',...
    [0.78 0.03 0.2 0.065],...
    'String',{'OD$_{87}$'},...
    'FontSize',fontsize-1,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'interpreter', 'latex',...
    'LineStyle','none',...
    'EdgeColor',[1 1 1]);


%%
set(gcf,'PaperPosition', [0 0 16 15]);
set(gcf,'PaperSize', [15.5 15]);

% print -dpdf buffer_Tstem_vs_OD87
% open buffer_Tstem_vs_OD87.pdf

%% Plotting: Optical density (OD) and T1 time

figure(16)
colormap(jet);
% subplot(2,1,2)
imagesc(h_vector,v_vector,OD); axis tight
% imagesc(h_vector,v_vector,T1_array/mus); axis tight
c_bar=colorbar;
set(c_bar,'FontSize',fontsize);
% caxis([0,  max(max(OD_array))])
caxis([0, 5]);

hold on
contourlines=linspace(0, 100,21);
[C h]=contour(h_vector,v_vector,spatial_res/mum,contourlines,'LineColor','w');
hold off
clabel(C,h,'Color','w')


set(gca,'YDir','normal','XDir','normal')
set(gca,'FontSize',fontsize);

xlabel('Temperature (^{\circ}C)');Hylabel=ylabel('Total P_{fill} (mbar)');
% str_title={sprintf('T1 (in \\mus) for various buffer pressures, with T=%g degC',T-273.15); 'Contour lines show the corresponding OD'};
% title('b) \deltaB_{photon}','FontSize',fontsize,'FontWeight','Normal');
title('Total OD, with spatial resolution contour lines','FontSize',fontsize,'FontWeight','Normal');

annotation('textbox',...
    [0.78 0.03 0.2 0.065],...
    'String',{'OD$_{87}$'},...
    'FontSize',fontsize-1,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'interpreter', 'latex',...
    'LineStyle','none',...
    'EdgeColor',[1 1 1]);
