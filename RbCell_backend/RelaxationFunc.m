function [T1, T2, gamma] = RelaxationFunc(cell_shape,cell_dimensions,T,Pbuffer_fill,Tfill,bgas,Rb)
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
%% buffer gas relaxation, diffusion constant
for i=1:length(Pbuffer_fill)
    t_bgas=bgas(i,:);
    
    v_av_bg(i)=sqrt(8*kB*T/pi/mass.(t_bgas));
    vrel_rbbg(i) = sqrt(Rb.v_av_rb^2+v_av_bg(i)^2);
    
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
I=3/2; % Rb87 nuclear spin
gamma.SE1 = Rb.gamma_rbrb_SE; %Characteristic time for spin exchange collisions
gamma.SE2 = (6*I+1)/(8*I+1)*gamma.SE1; %Characteristic time for spin exchange collisions

%% T1 and T2 times, combining all relaxation mechanisms
T1 = (gamma.SE1 + gamma.walls + gamma.bg1)^-1;
T2 = (gamma.SE2 + gamma.walls + gamma.bg2)^-1;
