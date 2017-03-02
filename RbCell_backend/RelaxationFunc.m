function [T1, T2, gamma] = RelaxationFunc(cell_shape,cell_dimensions,T,Pbuffer_fill,Tfill,bgas)
% Calculates T1 and T2 times, and the associated relaxation rates, for Rb87
% in an atomic vapor cell. The script can be adapted for other alkali
% species by inputting the appropriate constants.
% 
% gamma output: gamma.SE1, gamma.SE2, gamma.walls, gamma.bg1, gamma.bg2
% Rb is obtained from Rb=RbProperties(T,rb87_abundance)
% 
% see Andrew Horsley PhD thesis, Ch.2, for detailed theory: 
% High Resolution Field Imaging with Atomic Vapor Cells, Basel 2015
%
% Andrew Horsley 2017

constants; units; BufferGasProperties;

n_buffer= Pbuffer_fill/kB/Tfill; % derived buffer density, in m^3

%% Rb properties
I=3/2; % Rb87 nuclear spin
sigma_rbrb_SE = 1.9*10^-14*cm^2; %spin exchange cross section for rb, Walter (2002)

if T-273.15<39.3
    P= 10.^(2.881+4.857-4215./T)*torr; % Steck. Total Rb pressure with Rb in solid phase, in Pa
%     P = 133.323 * 10.^(-94.04826-1961.258./T-0.03771687*T+42.57526*log10(T));
else
    P = 10.^(2.881+4.312-4040./T)*torr; % Steck. Total Rb pressure with Rb in liquid phase, in Torr
%     P = 133.323 * 10.^(15.88253-4529.635./T+0.00058663*T-2.99138*log10(T));
end
density_total=P/kB./T; %Total Rb density at T, in m^3

v_av_rb= sqrt(8*kB.*T/pi/m87Rb); %Average Rb velocity
vrel_rbrb = sqrt(2)*v_av_rb; %Rb-Rb rel velocity
gamma_rbrb_SE=density_total*sigma_rbrb_SE.*vrel_rbrb; %SE collision rate

%% buffer gas relaxation, diffusion constant

for i=1:length(Pbuffer_fill)
    t_bgas=bgas(i,:);
    
    v_av_bg(i)=sqrt(8*kB*T/pi/mass.(t_bgas));
    vrel_rbbg(i) = sqrt(v_av_rb^2+v_av_bg(i)^2);
    
    t_gamma_bg1(i)=n_buffer(i)*sigma_1.(t_bgas)(1)*vrel_rbbg(i)*(T/sigma_1.(t_bgas)(2))^sigma_1.(t_bgas)(3); % buffer T1 relaxation rate
    t_gamma_bg2(i)=n_buffer(i)*sigma_2.(t_bgas)(1)*vrel_rbbg(i)*(T/sigma_2.(t_bgas)(2))^sigma_2.(t_bgas)(3); % buffer T2 relaxation rate
    t_D(i) =D0.(t_bgas)(1)*p_atmoshperic/ (Pbuffer_fill(i)*(D0.(t_bgas)(2)/Tfill)) *(T/D0.(t_bgas)(2))^D0.(t_bgas)(3); %Diffusion constant D, in m^2/s
    
end
gamma.bg1 = sum(t_gamma_bg1);
gamma.bg2 = sum(t_gamma_bg2);
diff_coeff = (sum(1./t_D))^-1;

%% Wall relaxation
% from diffusion equation in sphere, cylinder, rectangular box
% see Andrew Horsley PhD thesis, Ch.2

if strcmp(cell_shape ,'sphere')
    cell_radius = cell_dimensions;
    k2=(pi/cell_radius)^2;
elseif strcmp(cell_shape ,'cylinder')
    cell_radius = cell_dimensions(1);
    cell_thickness = cell_dimensions(2);
    mu1 = 2.405/cell_radius; %Using the actual cell size
    v1 = pi/cell_thickness;
    k2=(mu1^2+v1^2);
elseif strcmp(cell_shape ,'rectangle')
    a = cell_dimensions(1);
    b = cell_dimensions(2);
    c = cell_dimensions(3);
    k2=pi^2*(1/a^2+1/b^2+1/c^2);
end

gamma.walls = diff_coeff*k2;

%% Rb-Rb spin exchange relaxation (SE)
% From atomic freq. standards book (Vanier), the relationship
% between T1 and T2 spin-exchange relazation is given as 
% Gamma2 = (6I+1)/(8I+1) Gamma1 = 10/13 Gamma1 for Rb87
gamma.SE1 = gamma_rbrb_SE; %Characteristic time for spin exchange collisions
gamma.SE2 = (6*I+1)/(8*I+1)*gamma.SE1; %Characteristic time for spin exchange collisions

%% T1 and T2 times, combining all relaxation mechanisms
T1 = (gamma.SE1 + gamma.walls + gamma.bg1)^-1;
T2 = (gamma.SE2 + gamma.walls + gamma.bg2)^-1;
