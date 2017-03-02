function output = sensitivity_Bmw_func(T,Pfill,InputParam)
% Calculates the photon-shot-noise and atomic-projection-noise limited
% sensitivities for sensing microwave magnetic fields with an 87Rb vapor
% cell. 
% For the relevant theory, see:
% Horsley, Du and Treutlein, NJP 17 (2015) 112002; 
% Andrew Horsley PhD thesis, Ch.6:
% High Resolution Field Imaging with Atomic Vapor Cells, Basel 2015
%
% Andrew Horsley, 2017
constants; units;

Rb=RbProperties(T,InputParam.rb87_abundance);
[T1,T2,~]= RelaxationFunc(InputParam.cell_shape,InputParam.cell_dimensions,T,Pfill,InputParam.Tfill,InputParam.bgas);
BufferShiftBroadening = BufferShiftBroadeningFunc(T,Pfill,InputParam.Tfill,InputParam.bgas);

% mw pulse length that will be used for estimating the distance atoms 
% diffuse during a measurement:
dt_mw=T1; 
% The T1 time seems to give the best match with data, and is anyway
% essentially equal to T2 in thin cells (due to wall collisions). 

%% Calculating the optical depth and laser transmission:

if strcmp(InputParam.cell_shape ,'sphere')
    Loptical=InputParam.cell_dimensions(1);
elseif strcmp(InputParam.cell_shape ,'cylinder')
    Loptical=InputParam.cell_dimensions(2);
elseif strcmp(InputParam.cell_shape ,'rectangle')
    Loptical=InputParam.cell_dimensions(3);
end

% To account for the change in buffer gas shift as cell filling parameters 
% are varied, the buffer gas shift is automatically compensated for. The
% InputParam.Delta should be given assuming the buffer gas shift is zero
Delta = InputParam.detuning + BufferShiftBroadening.OptShift; 

n_buffer = Pfill/kB/InputParam.Tfill; % buffer density, in m^3
Esus_out=Esus_func(InputParam.Dline,InputParam.Bdc,T,Loptical,InputParam.rb87_abundance,Delta,n_buffer,InputParam.bgas);
OD=Esus_out.OD_p;
OD87=Esus_out.OD87_p;

I_transmitted=InputParam.I0*exp(-OD); % W/m2, laser intensity after cell

%% Diffusion-limited sensor (eg 50x50x100 um)
% The transverese spatial resolution limit is given by the distance atoms
% diffuse during a measurement (dtmw).

% diffusion-limited sensor size
diffusion_length = sqrt(2*dt_mw*BufferShiftBroadening.diff_coeff);
Vsensor=diffusion_length^2*Loptical; % diffusion limited sensor volume
Adiffusion = diffusion_length^2; % effective camera pixel size for a diffusion-limited sensor



%% Photon Shot Noise: Raw

% counts, shot noise per diffusion limited sensor
Nph_sensor = I_transmitted*Adiffusion*InputParam.Tprobe/hbar/om0; %# of photons / sensor / probe pulse
Nel_sensor = InputParam.QE*Nph_sensor; % number of electrons / sensor / probe pulse
delta_Nel_sensor = sqrt(Nel_sensor); % shot noise in the electron number
Ncounts = InputParam.gain*Nel_sensor; % registered camera counts
delta_Ncounts = InputParam.gain*delta_Nel_sensor; % noise in camera counts

%% Photon Shot Noise: ODmw noise
% We take actual and reference images, with/without the microwave present.
% The log of the actual/reference ratio is then the change in optical
% density induced by the microwave, ODmw:
% ODmw = -log(Nact/Nref) = 0 
% where Nact/Nref are the actual/reference image counts
% The actual and reference images counts are typically on the same order of
% magnitude. For simplicity in estimating the photon shot noise, we will
% assume that Nact=Nref=Ncounts. 
% The photon shot noise in ODmw is then
% delta_ODmw = sqrt( (dOD/dNact * deltaNact)^2 + (dOD/dNref * deltaNref)^2 )
%           = sqrt( 2*(dOD/Ncounts * delta_Ncounts)^2)

% partial derivative ODmw wrt Ncounts:
dODdNcounts = -1/Ncounts; 

% ODmw photon shot noise:
delta_ODmw = sqrt(2*(dODdNcounts*delta_Ncounts)^2);

% % ODmw photon shot noise, 
% delta_OD_shot_average=delta_ODmw/sqrt(Nshots);
% delta_OD_run_average=delta_ODmw/sqrt(Nshots)/sqrt(Nruns);

%% shot noise limited sensitivity
% Omega_min: smallest detectable Rabi frequency
% Bmw_min: corresponding smallest detectable mw magnetic field
% 
% deltaODmin = sqrt(hbar*om0/I0 * A * dt_im)
% Omega_min = dt_im/2 * sqrt(deltaODmin/OD)

% standard sequence (no pi/2 preparation)
Omega_min = 2/dt_mw*sqrt(delta_ODmw/(OD87*InputParam.TransPopDiff))*exp(dt_mw/T2);
Bmw_min = hbar/muB*Omega_min;
dBphoton = Bmw_min/sqrt(InputParam.DutyCycle);

% optimised sequence (preparation with pi/2 pulse)
Omega_min_opt = 2/dt_mw *delta_ODmw/(OD87*InputParam.TransPopDiff) * exp(dt_mw/T2);
Bmw_min_opt = hbar/muB*Omega_min_opt;
dBphoton_opt = Bmw_min_opt/sqrt(InputParam.DutyCycle);

%% Atomic Projection Noise

Tmeas = 1; % coherence time
atoms_interacting=1/3; % maximum fraction of prepared atoms that can interact with the mw (1/3 of F=1 atoms)
Nat = Rb.density_87*atoms_interacting*Vsensor;
dBatom = hbar/muB/sqrt(Nat*T2*Tmeas);

%% Output

output.T1=T1;
output.T2=T2;
output.delta_ODmw=delta_ODmw;
output.Omega_min=Omega_min;
output.Omega_min_opt=Omega_min_opt;
output.dBphoton=dBphoton;
output.dBphoton_opt=dBphoton_opt;
output.Nat=Nat;
output.dBatom=dBatom;
output.diffusion_length=diffusion_length;
output.Ncounts=Ncounts;
output.delta_Ncounts=delta_Ncounts;
output.OD=OD;
output.OD87=OD87;
end