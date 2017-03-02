% 4/7/14: Comparison of different Kr and N2 pressures in a buffer gas 
% mixtures, and the resulting T1 time
%
% 16/10/14 - fixed issue with N2 diffusion coefficient

clear
constants;
units;

%% Species of interest: 87Rb
rb87_abundance = 1; %natural abundance of 87Rb. Set to 1 for 100% 87Rb cell

%% Cell and laser beam dimensions
cell_thickness=200*mum; %thickness of rectangular cell (z-dimension)

% parameters
I0=30*mW/cm^2; % incident laser intensity


% temperature scan
T=linspace(80,170,10)+273.15; %cell temperature, in Kelvin
Tfill = 22+ 273.15; %standard temperature (eg buffer gas pressure is 10mbar @ T0) - 0degC from BOlsen thesis

% Buffer pressure scan
bgas=[sprintf('Kr'); sprintf('N2')]; %buffer gas types (as many as you want)
Ptotal = linspace(0,300,10)*mbar;
P1_frac = 0.75; P2_frac = 1-P1_frac;
P0_1 = P1_frac*Ptotal; P0_2 = P2_frac*Ptotal;
Pfill=[P0_1; P0_2];
n_buffer= Pfill/kB/Tfill; % buffer density, in m^3

%%   
%prepare arrays
T1_array=zeros(size(T,1),size(Pfill,1));
OD87_array=zeros(size(T,1),size(Pfill,1));
delta_array=zeros(size(T,1),size(Pfill,1));
Sensitivity=zeros(size(T,1),size(Pfill,1));
dBphoton=zeros(size(T,1),size(Pfill,1));
dBatom=zeros(size(T,1),size(Pfill,1));
spatial_res=zeros(size(T,1),size(Pfill,1));

for k1=1:size(T,2)
for k2=1:size(Ptotal,2)

    temp=T1_ultrathin_func(T(k1),rb87_abundance,cell_thickness,Pfill(:,k2),bgas);   
    T1_array(k2,k1) = temp.T1;
    delta_array(k2,k1) = temp.OptShift;
    
    temp2=sensitivity_ultrathin_func(T(k1),rb87_abundance,cell_thickness,Pfill(:,k2),bgas,delta_array(k2,k1),I0);
    Sensitivity(k2,k1)=temp2.Sensitivity;
    dBphoton(k2,k1)=temp2.Sensitivity_opt;
    dBatom(k2,k1)=temp2.dBatom;
    
    spatial_res(k2,k1) =temp2.diffusion_sensor_size;
    Omega_min(k2,k1)=temp2.Omega_min;
    Omega_min_opt(k2,k1)=temp2.Omega_min_opt;
    Nat(k2,k1)=temp2.Nat;
    delta_OD_run_average(k2,k1)=temp2.delta_OD_run_average;
    Ncounts(k2,k1)=temp2.Ncounts;
    delta_Ncounts(k2,k1)=temp2.delta_Ncounts;
    OD(k2,k1)=temp2.OD;
    OD87(k2,k1)=temp2.OD87;
end
end


%% Plotting
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

xlabel('Temperature (^{\circ}C)');Hylabel=ylabel('Total P_{fill} (mbar)');
title('a) \deltaB_{atom}','FontSize',fontsize,'FontWeight','Normal');

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

xlabel('Temperature (^{\circ}C)');Hylabel=ylabel('Total P_{fill} (mbar)');
title('b) \deltaB_{photon}','FontSize',fontsize,'FontWeight','Normal');

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

%%
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
title('OD_{87}','FontSize',fontsize,'FontWeight','Normal');

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

%%

figure(15)
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
% contour(v_vector,h_vector,OD_array,contourlines,'ShowText','on','LineColor','k')
[C h]=contour(h_vector,v_vector,spatial_res/mum,contourlines,'LineColor','w');
hold off
clabel(C,h,'Color','w')
% ,'Color','w'
% c_bar=colorbar;
% set(c_bar,'FontSize',fontsize);
% caxis([0,  max(max(OD_array))])

set(gca,'YDir','normal','XDir','normal')
set(gca,'FontSize',fontsize);

xlabel('Temperature (^{\circ}C)');Hylabel=ylabel('Total P_{fill} (mbar)');
% str_title={sprintf('T1 (in \\mus) for various buffer pressures, with T=%g degC',T-273.15); 'Contour lines show the corresponding OD'};
% title('b) \deltaB_{photon}','FontSize',fontsize,'FontWeight','Normal');
title('\deltaB_{photon}','FontSize',fontsize,'FontWeight','Normal');

annotation('textbox',...
    [0.78 0.03 0.2 0.065],...
    'String',{'OD$_{87}$'},...
    'FontSize',fontsize-1,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'interpreter', 'latex',...
    'LineStyle','none',...
    'EdgeColor',[1 1 1]);
