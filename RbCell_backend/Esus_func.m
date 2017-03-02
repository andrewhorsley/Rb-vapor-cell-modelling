function output=Esus_func(Dline,Bdc,T,L,purity,Delta,n_buffer,bgas)
% output=Esus_func(Dline,Bdc,T,L,purity,Delta,n_buffer,bgas)
% requires RbDlines_func
%
% Andrew Horsley, 25/1/2016
%
% Electric susceptibility for light travelling through a vapor cell
% The (complex) susceptibility, Esus, is calculated separately for sigma+/- 
% and pi polarised light. The real component of Esus gives the dispersion,
% whilst the imaginary part gives the absorption
% I.e., the output is slightly different to true electrical susceptibility:
% Here, the optical density is OD=imag(Esus) and dispersion=real(Esus), 
% instead of the usual OD=kL*imag(Esus) and dispersion=kL*real(Esus), where
% kL is the laser wavevector
%
% Based off the Durham model from the Adams/Hughes group, initially
% presented in [1] and extended by L. Weller [2] and M. Zentile [3]. 
% The output has been checked against the python program Elecsus [4],
% which in turn has shown close agreement with experiment.
%
% The model assumes the laser beam is weak, and therefore neglects optical
% pumping and saturation effects. The model includes:
% Buffer gas collisional broadening and shifts, for any mixture of gases
% listed in BufferGasProperties.m
% Rb self-broadening, due to dipole-dipole interactions at high tempertures
% (and thus high Rb densities) [5]
% The effect of arbitrarily large dc magnetic fields (up to e.g. several
% Tesla, into the hyperfine Paschen-Back regime) [6]
% Note that the model breaks for Bdc=0 (ElecSus has the same problem). We
% get around this by calculating using Bdc=1e-6 T when Bdc=0 is input.
%
% INPUTS:
% Dline is either 'D1' or 'D2'. 
% Bdc is the dc magnetic field [T], which can be a vector
% T is the cell temperature [K], L the cell length [m], purity the fraction
% of 87Rb atoms in the cell (0.2783 = natural abundance)
% Delta is the detuning (angular frequency) with respect to the weighted D1
% or D2 line centre (zero value in Steck). Delta can also be input as a vector
% n_buffer is a vector containing the densities of each buffer gas in the
% cell, and bgas is a string vector listing the corresponding buffer gas
% species. He, Ne, N2, Ar, Kr, and Xe are currently supported (25/1/2016)
%
% Checked against ElecSus at Bdc=0: good +/-0.5%
% larger deviations at larger Bdc
%
% [1] P. Siddons et al., (2008) J. Phys. B 41(55) 155004
% [2] L. Weller, PhD Thesis, (2013), Durham University
% [3] M. Zentile, PhD Thesis (2015), Durham University
% [4] M. Zentile et al., (2015) ElecSus: A program to calculate the electric
% susceptibility of an atomic ensemble. Computer Physics Communications, 189, 162–174. 
% [5] L. Weller er al., (2011) J. Phys. B 44(19) 195006
% [6] L. Weller er al., (2012) J. Phys. B 45(21) 215005 

BufferGasProperties; 

kB=1.3806504e-23; % Boltzmann's Constant
hbar=1.054571628e-34; % Planck's Constant/2pi
% echarge=1.602176487e-19; % Elementary Charge, pre 26/1/2016
echarge=1.602176565e-19; % Elementary Charge, ElecSus value
eps0=8.854187817e-12; % Permittivity of Vacuum
% a0=0.52917720859e-10; % Bohr Radius, pre 26/1/2016
a0=0.52917721092e-10; % Bohr Radius, ElecSus value
m87Rb=1.443160648e-25; % Rb87 Atomic Mass
m85Rb=1.409993199e-25; % Rb85 Atomic Mass
torr=133.3224;
cm=1e-2;
GHz=1e9;

I87=3/2; % 87Rb nuclear spin
I85=5/2; % 85Rb nuclear spin
J_5S1_2 = 1/2; % 5S1/2 electron
J_5P1_2 = 1/2; % 5P1/2 electron
J_5P3_2 = 3/2; % 5P3/2 electron
if strcmp(Dline,'D1')
%     lambda=794.978851156e-9; % Wavelength D2 Transition (Vacuum, Steck)
    lambda=794.978969380e-9; % Wavelength D2 Transition (Vacuum, ElecSus)
    dred=5.182*echarge*a0; % D1 reduced matrix element <Lg=0||er||Le=1>
    Jfactor = (2*J_5S1_2+1)/(2*J_5P1_2+1); % normalisation factor of (2J+1)/(2J'+1) = 1 for D1, 1/2 for D2
    Gamma=2*pi*5.746e6; %Rb85,87 Natural Line Width D1 Transition (ElecSus)
    BetaDipole = 2*pi*0.69*10^-7*cm^3; % Rb D2 dpiole-dipole broadening, L. Weller, et al, J. Phys. B 45, 215005 (2012) 
elseif strcmp(Dline,'D2')
%     lambda=780.241209686e-9; % Wavelength D2 Transition (Vacuum, Steck)
    lambda=780.24132411e-9; % Wavelength D2 Transition (Vacuum, ElecSus)
    dred=5.177*echarge*a0; % D2 reduced matrix element <Lg=0||er||Le=1>
    Jfactor = (2*J_5S1_2+1)/(2*J_5P3_2+1); % normalisation factor of (2J+1)/(2J'+1) = 1 for D1, 1/2 for D2
    Gamma=2*pi*6.0666e6; % Rb85,87 Natural FWHM Line Width D2 Transition (Steck)    
%     Gamma=2*pi*6.065e6; % Rb85,87 Natural FWHM Line Width D2 Transition (ElecSus)   
    BetaDipole = 2*pi*1.10*10^-7*cm^3; % Rb D2 dpiole-dipole broadening, L. Weller, et al, J. Phys. B 45, 215005 (2012)
end
u87=sqrt(2*kB*T/m87Rb); % 1/e width of Maxwell-Boltzmann-Gaussian in 1D
u85=sqrt(2*kB*T/m85Rb); % 1/e width of Maxwell-Boltzmann-Gaussian in 1D
kL=2*pi/lambda; % Wave Number D2 Transition (Vacuum)

%% Rb vapor density

% Rb vapor pressure (Steck, alkali data)
pRb = zeros(size(T));
for nn=1:length(T)
    if T(nn)<273.15+39.30 % vapor pressure solid phase (Pa)
%     pRb(nn) = 133.32 * 10.^(2.881+4.857-4215./T(nn)); % Steck
    pRb(nn) = 133.32 * 10.^(2.881+4.857-4215./T(nn)); % ElecSus
    else % vapor pressure liquid phase (Pa)
%     pRb(nn) = 133.32 * 10.^(2.881+4.312-4040./T(nn)); % Steck
    pRb(nn) = 133.32 * 10.^(2.881+8.316-4275./T(nn)-1.3102*log10(T(nn))); % ElecSus
    end
end

nRb = pRb./(kB*T); % density of Rb atoms
nRb87 = purity*nRb; % density of 87Rb (for natural abundance 27.83%)
nRb85 = (1-purity)*nRb; % density of 85Rb

%% Optical broadening and lineshift
% See Ch.2 of PhD thesis, Andrew Horsley, University of Basel, 2015:
% High Resolution Field Imaging with Atomic Vapor Cells
%%%%%% DIFFERENT D1 and D2 collisional broadening, lineshifts???
% Currently using D2 parameters for the D1 broadening and lineshift

T0degC=273.15;
P_buffer0 = n_buffer*kB*T0degC;
for i=1:length(n_buffer)
    t_bgas=bgas(i,:);
    
    t_OptBroad(i,:)=opt_broad.(t_bgas)(1)*P_buffer0(i).*(T/opt_broad.(t_bgas)(2)).^opt_broad.(t_bgas)(3)*(opt_broad.(t_bgas)(2)/T0degC); %optical buffer broadening
    t_OptShift(i,:)=opt_shift.(t_bgas)(1)*P_buffer0(i).*(T/opt_shift.(t_bgas)(2)).^opt_shift.(t_bgas)(3)*(opt_shift.(t_bgas)(2)/T0degC); %optical line shift
end

% combine properties of the (one or many) buffer gases as single quantities:
Gammapress = sum(2*pi*t_OptBroad); % Rb D2 broadening from buffer gas
deltapress = sum(2*pi*t_OptShift); % Rb D2 line shift from buffer gas

% Rb D2 dpiole-dipole broadening:
% Gamma_dipole = Beta*N; Beta1 = 2pi*0.69*10^-7 Hz cm3, Beta2 = 2pi*1.10*10^-7 Hz cm3, N is Rb density
% L. Weller, et al, J. Phys. B 45, 215005 (2012) 
Gamma_dipole = BetaDipole*nRb;  

% add natural, pressure-broadened, and dipole-dipole linewidths (all Lorentzian)
Gammatot = Gamma + Gammapress + Gamma_dipole;
a87=Gammatot./(kL.*u87);
a85=Gammatot./(kL.*u85);

%% Transition strengths and frequencies
% The RbDlines function calculates the energy levels, optical transition
% strenghts and optical transition frequencies
% RbDlines_func can accept Bdc as a vector
RbDlines=RbDlines_func(Dline,Bdc);

% frequency offset of 85Rb Fg->Fe transitions
dom85p=2*pi*RbDlines.TransFreq85_p + deltapress; % pi polarised transitions
dom85m=2*pi*RbDlines.TransFreq85_m + deltapress; % sigma- polarised transitions
dom85z=2*pi*RbDlines.TransFreq85_z + deltapress; % sigma+ polarised transitions

% frequency offset of 87Rb Fg->Fe transitions
dom87p=2*pi*RbDlines.TransFreq87_p+ deltapress; 
dom87m=2*pi*RbDlines.TransFreq87_m+ deltapress; 
dom87z=2*pi*RbDlines.TransFreq87_z+ deltapress; 

% 85Rb Transition strengths
CF85p = Jfactor*(RbDlines.TransStrength85_p).^2;
CF85m = Jfactor*(RbDlines.TransStrength85_m).^2;
CF85z = Jfactor*(RbDlines.TransStrength85_z).^2;
% 87Rb Transition strengths
CF87p = Jfactor*(RbDlines.TransStrength87_p).^2;
CF87m = Jfactor*(RbDlines.TransStrength87_m).^2;
CF87z = Jfactor*(RbDlines.TransStrength87_z).^2;

%% lineshape function for 87Rb and 85Rb
% Electrical susceptibility (chi) of the atoms for a given transition
% real part leads to dispersion index, imaginary part leads to OD
% It's a normalized Voigt profile centered at y=0 for 87Rb and 85Rb
% Based off Siddons 2008

sI87 = @(y)(sqrt(pi)*real(Faddeeva_erfcx(a87/2-1i*y)));
sR87 = @(y)(-sqrt(pi)*imag(Faddeeva_erfcx(a87/2-1i*y)));
chiF87 = @(y,CF2,dred)(CF2*dred^2*nRb87/(2*(2*I87+1))/(hbar*eps0)*(sR87(y)+1i*sI87(y))/(kL*u87));

sI85 = @(y)(sqrt(pi)*real(Faddeeva_erfcx(a85/2-1i*y)));
sR85 = @(y)(-sqrt(pi)*imag(Faddeeva_erfcx(a85/2-1i*y)));
chiF85 = @(y,CF2,dred)(CF2*dred^2*nRb85/(2*(2*I85+1))/(hbar*eps0)*(sR85(y)+1i*sI85(y))/(kL*u85));

%% sum over transitions with corresponding transition strengths (CF2) and detunings

% initiate output arrays:
output.Disp_p=zeros(length(Delta),length(Bdc)); 
output.Disp_m=zeros(length(Delta),length(Bdc)); 
output.Disp_z=zeros(length(Delta),length(Bdc)); 
output.OD_p=zeros(length(Delta),length(Bdc)); 
output.OD_m=zeros(length(Delta),length(Bdc)); 
output.OD_z=zeros(length(Delta),length(Bdc)); 
output.Disp87_p=zeros(length(Delta),length(Bdc)); 
output.Disp87_m=zeros(length(Delta),length(Bdc)); 
output.Disp87_z=zeros(length(Delta),length(Bdc)); 
output.OD87_p=zeros(length(Delta),length(Bdc)); 
output.OD87_m=zeros(length(Delta),length(Bdc)); 
output.OD87_z=zeros(length(Delta),length(Bdc)); 
output.Disp85_p=zeros(length(Delta),length(Bdc)); 
output.Disp85_m=zeros(length(Delta),length(Bdc)); 
output.Disp85_z=zeros(length(Delta),length(Bdc)); 
output.OD85_p=zeros(length(Delta),length(Bdc)); 
output.OD85_m=zeros(length(Delta),length(Bdc)); 
output.OD85_z=zeros(length(Delta),length(Bdc)); 

for bb=1:length(Bdc)
chi87p=zeros(1,length(Delta)); chi87m=zeros(1,length(Delta)); chi87z=zeros(1,length(Delta));
chi85p=zeros(1,length(Delta)); chi85m=zeros(1,length(Delta)); chi85z=zeros(1,length(Delta));
for k=1:size(dom87p,1)
    chi87p = chi87p + chiF87( (Delta-dom87p(k,bb))./(kL.*u87) , CF87p(k,bb),dred );
    chi87m = chi87m + chiF87( (Delta-dom87m(k,bb))./(kL.*u87) , CF87m(k,bb),dred );
end
for k=1:size(dom87z,1)
    chi87z = chi87z + chiF87( (Delta-dom87z(k,bb))./(kL.*u87) , CF87z(k,bb),dred );
end
for k=1:size(dom85p,1)
    chi85p = chi85p + chiF85( (Delta-dom85p(k,bb))./(kL.*u85) , CF85p(k,bb),dred );
    chi85m = chi85m + chiF85( (Delta-dom85m(k,bb))./(kL.*u85) , CF85m(k,bb),dred );
end
for k=1:size(dom85z,1)
    chi85z = chi85z + chiF85( (Delta-dom85z(k,bb))./(kL.*u85) , CF85z(k,bb),dred );
end

% Phase shift due to vapor cell (.: includes L)
output.Disp_p(:,bb) = 1/2*real(chi87p+chi85p)*kL*L;
output.Disp_m(:,bb) = 1/2*real(chi87m+chi85m)*kL*L;
output.Disp_z(:,bb) = 1/2*real(chi87z+chi85z)*kL*L;

% Optical density of vapor cell (.: includes L)
output.OD_p(:,bb) = imag(chi87p+chi85p)*kL*L;
output.OD_m(:,bb) = imag(chi87m+chi85m)*kL*L;
output.OD_z(:,bb) = imag(chi87z+chi85z)*kL*L;

% Same for individual 85Rb and 87Rb contributions
output.Disp87_p(:,bb) = 1/2*real(chi87p)*kL*L;
output.Disp87_m(:,bb) = 1/2*real(chi87m)*kL*L;
output.Disp87_z(:,bb) = 1/2*real(chi87z)*kL*L;
output.OD87_p(:,bb) = imag(chi87p)*kL*L;
output.OD87_m(:,bb) = imag(chi87m)*kL*L;
output.OD87_z(:,bb) = imag(chi87z)*kL*L;

output.Disp85_p(:,bb) = 1/2*real(chi85p)*kL*L;
output.Disp85_m(:,bb) = 1/2*real(chi85m)*kL*L;
output.Disp85_z(:,bb) = 1/2*real(chi85z)*kL*L;
output.OD85_p(:,bb) = imag(chi85p)*kL*L;
output.OD85_m(:,bb) = imag(chi85m)*kL*L;
output.OD85_z(:,bb) = imag(chi85z)*kL*L;
end

