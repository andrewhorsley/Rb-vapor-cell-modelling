% properties of common buffer gases and their interaction with Rb. 
% see Ch.2 of my thesis for a description of the theory
% https://atom.physik.unibas.ch/publications/pdf/PhD_Thesis_Horsley_Andrew.pdf
%
% Andrew Horsley 27.5.2015

% units;
cm=1e-2;
MHz=1e6;
torr=133.3224;

% mass (1 amu = 1.66053892*10^-27 kg)
mass.He=6.6465e-27; %mass natHe, using mass natHe=4.002602 amu (google)
mass.Ne=3.3509*10^-26; %mass natNe, using mass natNe=20.1797 amu (google)
mass.N2=2*2.3259*10^-26; %mass natN2, using mass natN=14.0067 amu (google)
mass.Ar=6.6335e-26; %mass natAr, using mass=39.948 amu (google)
mass.Kr=1.3915e-25; %mass natKr, using mass=83.798 amu (google)
mass.Xe=2.1802e-25; %mass natXe, using mass=131.293 amu (google)

% sigma total: sigma, T0, kappa
sigma_total.He = [1.52e-14*cm^2/20,   0+273.15,  1/2]; % He Croucher 1969
sigma_total.Ne = [268*10^-16*cm^2/20, 0+273.15,  1/2]; % Ne Croucher 1969
sigma_total.N2 = [785*10^-16*cm^2/20, 0+273.15,  1/2]; % N2 Croucher 1969
sigma_total.Ar = [5.72e-14*cm^2/20,   0+273.15   1/2]; % Ar Croucher 1969
sigma_total.Kr = [1e-13*cm^2/20,      0+273.15   1/2]; % Kr Croucher 1969
sigma_total.Xe = [0,          0+273.15   1/2]; % Xe - no data available

% sigma T1: sigma, T0, kappa
sigma_1.He = [8.7e-24*cm^2,    150+273.15, 0.42]; % He Happer2010
sigma_1.Ne = [1.9e-23*cm^2,    30+273.15,  0.27]; % Ne Vanier1989
sigma_1.N2 = [1e-22*cm^2,      70+273.15,  0.3]; % N2 Happer2010
sigma_1.Ar = [6.1e-22*cm^2,    27+273.15   0.32]; % Ar Happer2010
sigma_1.Kr = [2.7e-21*cm^2,    27+273.15   0.32]; % Kr Happer2010
sigma_1.Xe = [2.2e-19*cm^2,    150+273.15  0.32]; % Xe Happer2010

% sigma T2: sigma, T0, kappa
sigma_2.He = [2.94e-21*cm^2,    27+273.15, 0.42]; % He Vanier1974
sigma_2.Ne = [5.55e-21*cm^2,    27+273.15,  0.27]; % Ne Vanier1974
sigma_2.N2 = [7.43e-21*cm^2,    27+273.15,  0.3]; % N2 Vanier1974 
sigma_2.Ar = [3.71e-21*cm^2,    27+273.15   0.32]; % Ar Vanier1974
sigma_2.Kr = [2.7e-20*cm^2,     27+273.15   0.32]; % Kr guess
sigma_2.Xe = [1.1e-18*cm^2,     21+273.15  0.32]; % Xe Brattke1998

% diffusion constant: D0, T0, kappa
D0.He = [0.35*cm^2  80+273.15  3/2]; % Nelson2001: 0.238cm2/s at 0degC
D0.Ne = [0.2*cm^2  0+273.15  3/2]; % Chrapkiewicz2014 - agree with Parniak2014
D0.N2 = [0.159*cm^2  60+273.15  3/2]; % Ishikawa2000 0.118cm2/s at 0degC
D0.Ar = [0.14*cm^2  27+273.15  3/2]; % Vanier1974 0.1215cm2/s at 0degC
D0.Kr = [0.068*cm^2  0+273.15  3/2]; % Chrapkiewicz2014 - agree with Parniak2014
D0.Xe = [0.057*cm^2  0+273.15  3/2]; % Chrapkiewicz2014 - agree with Parniak2014

% optical broadening: delta, T0, kappa
% kappa values are theoretical, from Kielkopf1976
opt_broad.He = [0*MHz/torr  0  0.42];
opt_broad.Ne = [9.47*MHz/torr  394  0.27]; % Rotondaro1997
opt_broad.N2 = [18.3*MHz/torr  394  0.3]; % Rotondaro1997
opt_broad.Ar = [17.7*MHz/torr  394  0.32]; % Rotondaro1997
opt_broad.Kr = [17.2*MHz/torr  394  0.32]; % Rotondaro1997
opt_broad.Xe = [17.8*MHz/torr  394  0.32]; % Rotondaro1997

% optical shift: delta, T0, kappa
% kappa values are theoretical, from Kielkopf1976
opt_shift.He = [0*MHz/torr  0  1.20];
opt_shift.Ne = [-2.44*MHz/torr  394  0]; % Rotondaro1997
opt_shift.N2 = [-5.79*MHz/torr  394  0.3]; % Rotondaro1997
opt_shift.Ar = [-5.76*MHz/torr  394  0.31]; % Rotondaro1997
opt_shift.Kr = [-5.50*MHz/torr  394  0.31]; % Rotondaro1997
opt_shift.Xe = [-6.19*MHz/torr  394  0.31]; % Rotondaro1997

% hyperfine shift: beta, delta, gamma, T0
mw_shift.He = [0/torr 0/torr 0/torr 333];
mw_shift.Ne = [392/torr  0/torr 0/torr 333];
mw_shift.N2 = [546.9/torr 0.55/torr 1.5e-3/torr 333];
mw_shift.Ar = [-59.7/torr -0.32/torr -3.5e-3/torr 333];
mw_shift.Kr = [-593.5/torr -0.47/torr 0/torr 333];
mw_shift.Xe = [0/torr 0/torr 0/torr 333];

% quenching cross section: D1, D2, T0, kappa
sigmaQ.He = [0 0 0 0];
sigmaQ.Ne = [0 0 0 0];
sigmaQ.N2 = [5.8e-15*cm^2, 4.3e-15*cm^2, 273.15, 0];
sigmaQ.Ar = [0 0 0 0];
sigmaQ.Kr = [0 0 0 0];
sigmaQ.Xe = [0 0 0 0];