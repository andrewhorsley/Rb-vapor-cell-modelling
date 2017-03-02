function OD=ODfuncPress(T,L,purity,Delta,n_buffer,bgas)
% OD = ODfuncPress(T,L,purity,Delta,n_buffer,buffergas)
%
% optical density of Rb in vapor cell assuming no pumping and no
% saturation, but including collisional broadening in Ne buffer gas
% see P. Siddons et al., J. Phys. B: At. Mol. Opt. Phys. 41, 155004 (2008)
% modified by including pressure broadening as in 
% M. D. Rotondaro et al., J. Quant. Spectrosc. Radiat. Transfer 57, 497 (1997)
%
% T is the cell temperature [K], L the cell length [m], purity the fraction
% of 87Rb atoms in the cell (0.2783 = natural abundance), and Delta is the
% detuning (angular frequency) with respect to the 87Rb Fg=2->Fe=2 line center.
% The Ne buffer gas pressure is pBuffer [Pa].
% 
% All hyperfine components of D2 line of 87Rb and 85Rb are considered
%
% AJH 28.5.2013 - modified to include Rb dipole-dipole interactions, as in
% L. Weller, et al, Journal of Physics B: Atomic, Molecular and Optical Physics 45, 215005 (2012) 
% Gamma_dipole = Beta*N; Beta1 = 2pi*0.69*10^-7 Hz cm3, Beta2 = 2pi*1.10*10^-7 Hz cm3, N is Rb density
%
% Update 27/2/2014: Changed cell_data and ODfuncPress so that they are now
% compatible with using one or multiple buffer gas types at once
% The frequency offsets (dom85, dom87) now include the lineshift due to the
% buffer gas/es
%
% Update 29/6/2015: Noted that erfc = 1-erf breaks for large detunings.
% Replaced erfz with Faddeeva package, in general scripts
% Still haven't worked out why the script breaks for large detunings
%
% Update 25/1/2016: Changed the exp(z^2)erfc(z) terms in sI(y) to erfcx(z),
% which avoids the problem of the solution disappearing for large
% detunings
% ODfuncPress can now also accept Delta as a vector input

kB=1.3806504e-23; % Boltzmann's Constant
hbar=1.054571628e-34; % Planck's Constant/2pi
echarge=1.602176487e-19; % Elementary Charge
eps0=8.854187817e-12; % Permittivity of Vacuum
a0=0.52917720859e-10; % Bohr Radius
m87Rb=1.443160648e-25; % Rb87 Atomic Mass
m85Rb=1.409993199e-25; % Rb85 Atomic Mass
Gamma=2*pi*6.0666e6; % Rb87 Natural Line Width D2 Transition (FWHM)
lambda=780.241209686e-9; % Wavelength D2 Transition (Vacuum)
kL=2*pi/lambda; % Wave Number D2 Transition (Vacuum)
torr=133.3224;
BufferGasProperties;
%% parameters of the cell (all SI units)

% Rb vapor pressure (Steck, alkali data *or* Siddons et al., J. Phys. B: AMOP 41, 155004 (2008))
pRb = zeros(size(T));
for nn=1:length(T)
    if T(nn)<273.15+39.30 % vapor pressure solid phase (Pa)
    pRb(nn) = 133.32 * 10.^(2.881+4.857-4215./T(nn)); % Steck
%     pRb(nn) = 133.323 * 10.^(-94.04826-1961.258./T(nn)-0.03771687*T(nn)+42.57526*log10(T(nn)));
    else % vapor pressure liquid phase (Pa)
    pRb(nn) = 133.32 * 10.^(2.881+4.312-4040./T(nn)); % Steck
%     pRb(nn) = 133.323 * 10.^(15.88253-4529.635./T(nn)+0.00058663*T(nn)-2.99138*log10(T(nn)));
    end
end

nRb = pRb./(kB*T); % density of Rb atoms
nRb87 = purity*nRb; % density of 87Rb (for natural abundance 27.83%)
nRb85 = (1-purity)*nRb; % density of 85Rb

%% absorption coefficient for 87Rb and 85Rb
I87=3/2; % 87Rb nuclear spin
I85=5/2; % 85Rb nuclear spin
dred=5.177*echarge*a0; % reduced matrix element <Lg=0||er||Le=1>
u87=sqrt(2*kB*T/m87Rb); % 1/e width of Maxwell-Boltzmann-Gaussian in 1D
u85=sqrt(2*kB*T/m85Rb); % 1/e width of Maxwell-Boltzmann-Gaussian in 1D

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

% Gamma_dipole = Beta*N; Beta1 = 2pi*0.69*10^-7 Hz cm3, Beta2 = 2pi*1.10*10^-7 Hz cm3, N is Rb density
% L. Weller, et al, Journal of Physics B: Atomic, Molecular and Optical Physics 45, 215005 (2012) 
cm=1e-2;
Beta2 = 2*pi*1.10*10^-7*cm^3;
Gamma_dipole = Beta2*nRb;  % Rb D2 dpiole-dipole broadening

% Gammatot = Gamma + Gammapress;
Gammatot = Gamma + Gammapress + Gamma_dipole;% add natural, pressure-broadened, and dipole-dipole linewidths (Lorentzian)
a87=Gammatot./(kL.*u87);
a85=Gammatot./(kL.*u85);

%Delta = 2*pi*[-2:0.01:2]*1e9; % detuning vector (Delta = 0 corresponds to 87Rb Fg=2->Fe=2 line center)

% frequency offset of 87Rb Fg->Fe transitions, order: 21,22,23,10,11,12
dom87 = 2*pi*[-2735.05; -2578.11; -2311.26; 4027.403; 4099.625; 4256.57]*1e6 + deltapress; 
% frequency offset of 85Rb Fg->Fe transitions, order: 32,33,34,21,22,23
dom85 = 2*pi*[-1371.29; -1307.87; -1186.91; 1635.454; 1664.714; 1728.134]*1e6 + deltapress; 

% HF transition strength CF2 for linear pol., order: 21,22,23,10,11,12
% same strength for other light polarizations as long as all atomic mF levels equally populated
CF87 = [1/18; 5/18; 7/9; 1/9; 5/18; 5/18];
% HF transition strength CF2 for linear pol., order: 32,33,34,21,22,23
% same strength for other light polarizations as long as all atomic mF levels equally populated
CF85 = [10/81; 35/81; 1; 1/3; 35/81; 28/81];

% normalized Voigt profile centered at y=0 for 87Rb and 85Rb
% sI87 = @(y) ( sqrt(pi)*real( exp(1/4*(a87-1i*2*y).^2).*(Faddeeva_erfc(a87/2-1i*y)) ) ); % pre 160125 - blows up for large detunings
% sI85 = @(y) ( sqrt(pi)*real( exp(1/4*(a85-1i*2*y).^2).*(Faddeeva_erfc(a85/2-1i*y)) ) ); % pre 160125 - blows up for large detunings
sI87 = @(y)(sqrt(pi)*real(Faddeeva_erfcx(a87/2-1i*y)));
sR87 = @(y)(-sqrt(pi)*imag(Faddeeva_erfcx(a87/2-1i*y)));

sI85 = @(y)(sqrt(pi)*real(Faddeeva_erfcx(a85/2-1i*y)));
sR85 = @(y)(-sqrt(pi)*imag(Faddeeva_erfcx(a85/2-1i*y)));

% absorption coeff for given transition Fg->Fe, linear polarization,
% detuning y=Delta/(kL*u) from line center, 
% strength (sqared matrix el.) CF2 (in units of reduced matrix el. dred^2).
alphaF87 = @(y,CF2) ( kL*CF2*dred^2.*nRb87.*1/(2*(2*I87+1)).*1/(hbar*eps0).*sI87(y)./(kL*u87) );
alphaF85 = @(y,CF2) ( kL*CF2*dred^2.*nRb85.*1/(2*(2*I85+1)).*1/(hbar*eps0).*sI85(y)./(kL*u85) );

% sum over transitions with corresponding CF2 and detunings
alpha87 = zeros(1,length(Delta));
alpha85 = zeros(1,length(Delta));
for k=1:length(dom87)
    alpha87 = alpha87 + alphaF87( (Delta-dom87(k)-2*pi*2578.11e6)./(kL.*u87) , CF87(k) );
    alpha85 = alpha85 + alphaF85( (Delta-dom85(k)-2*pi*2578.11e6)./(kL.*u85) , CF85(k) );
end

% optical density
OD = alpha87*L + alpha85*L; 
% OD = alpha87*L




