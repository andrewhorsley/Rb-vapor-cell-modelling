% fundamental constants in SI-units (taken from D.A. Steck, 
% "Rb 87 D Line Data", revision 2.1.4 (2010), Los Alamos)
%
% Andrew Horsley

% Fundamental Physical Constants (CODATA 1998)

clight=2.99792458e8; % Speed of Light
mu0=4*pi*1e-7; % Permeability of Vacuum
eps0=8.854187817e-12; % Permittivity of Vacuum
hplanck=6.62606896e-34; % Planck's Constant
hbar=1.054571628e-34; % Planck's Constant/2pi
echarge=1.602176487e-19; % Elementary Charge
muB=9.27400915e-24; % Bohr Magneton
me=9.10938215e-31; % Electron Mass
a0=0.52917720859e-10; % Bohr Radius
kB=1.3806504e-23; % Boltzmann's Constant
eV=1.60217657e-19; %eV in joules
Ry = 13.60569253*eV; % Rydberg constant, in joules

% 87Rb Physical Properties

m87Rb=1.443160648e-25; % Atomic Mass
om0=2*pi*384.2304844685e12; % Frequency D2 Transition
lambda=780.241209686e-9; % Wavelength D2 Transition (Vacuum)
kL=2*pi/lambda; % Wave Number D2 Transition (Vacuum)
vrec=hbar*kL/m87Rb; % Recoil Velocity D2
Erec=(hbar*kL)^2/(2*m87Rb); % Recoil Energy D2
omrec=hbar*kL^2/(2*m87Rb); % Recoil Frequency D2
Trec=2*Erec/kB; % Recoil Temperature D2
Gamma=2*pi*6.0666e6; % Natural Line Width D2 Transition (FWHM)
Isat=pi*hplanck*clight*Gamma/(3*lambda^3); % Saturation Intensity for Cooling Transition
as=100.44*a0; % S-wave scattering length for state |F=1,m_F=-1>, 
% for |2,1>: as=95.47*a0, between the states: as=98.09*a0, see Harber et al., PRA 66, 053616 (2002).
alpha0=hplanck*0.0794e-4; % Ground state DC electric polarizability
Ahfs = hplanck*3.41734130545215e9; % 87Rb 5S1/2 hyperfine constant;
Ehfs = hplanck*6.834682610904312e9; % 87Rb 5S1/2 hyperfine splitting;
