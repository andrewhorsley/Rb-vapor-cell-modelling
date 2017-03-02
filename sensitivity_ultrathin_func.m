function output = sensitivity_ultrathin_func(T,rb87_abundance,cell_thickness,Pfill,bgas,Delta,I0)
constants; units; BufferGasProperties;

Tfill=22+273.15;
n_buffer= Pfill/kB/Tfill; % buffer density, in m^3

temp=T1_ultrathin_func(T,rb87_abundance,cell_thickness,Pfill,bgas);
T1=temp.T1;
diff_coeff=temp.diff_coeff;
% diffusion_sensor_size = sqrt(2*T1*diff_coeff);
diffusion_sensor_size = 42*mum;

%% Calculating Rb properties:
% rb87_abundance = 0.2783; %natural abundance of 87Rb. Set to 1 for 100% 87Rb cell
if T-273.15<39.3
    P_temp = 10^(2.881+4.857-4215/T); % Steck. Total Rb pressure with Rb in solid phase, in Torr
%     P_temp = 133.323 * 10.^(-94.04826-1961.258./T-0.03771687*T+42.57526*log10(T));
    P=P_temp*torr; %pressure in Pa
else
    P_temp = 10^(2.881+4.312-4040/T); % Steck. Total Rb pressure with Rb in liquid phase, in Torr
%     P_temp = 133.323 * 10.^(15.88253-4529.635./T+0.00058663*T-2.99138*log10(T));
    P=P_temp*torr; %pressure in Pa
end

density_totRb=P/kB/T; %Total Rb density at T, in m^3
density_87Rb=density_totRb*rb87_abundance; %87Rb density at T, in m^3

% temp2 = RbProperties(T,rb87_abundance);
% density_totRb=temp2.density_total;
% density_87Rb=temp2.density_87;

%% Calculating the optical depth:
OD=ODfuncPress(T,cell_thickness,rb87_abundance,Delta,n_buffer,bgas);
OD87=ODfuncPress(T,cell_thickness,1,Delta,n_buffer,bgas)*rb87_abundance; %OD due to Rb87 only

%% Photon Shot Noise
% per diffusion-limited sensor (eg 50x50x100 um)

I_transmitted=I0*exp(-OD); % W/m2, laser intensity after cell
Tprobe = 0.3*mus; %s, probe pulse length
% diffusion_sensor_size = 42*mum; %size of a diffusion limited sensor
Vsensor=diffusion_sensor_size^2*cell_thickness; % diffusion limited sensor volume
% Apixel = (5.6*mum)^2; %m2, pixel area
% M = 0.28; %image magnification - 20um at cell to 5.6um at ccd
% Aeff = Apixel/M^2; % effective pixel size
Adiffusion = diffusion_sensor_size^2; % effective pixel size for a diffusion-limited sensor
QE = 0.27; % quantum efficiency of Guppy Pro camera
gain=1; % default set to 0dBm in matcam_lite
optical_pumping=1; %optical pumping efficiency

% fprintf('set: %g um cell \n',cell_thickness/mum);
% fprintf('set: Tstem = %g degC \n',T-273.15);
% fprintf('OD = %g \n',OD);
% fprintf('OD-87 = %g \n',OD87);
% fprintf('ODmw = %g \n',OD87*optical_pumping/3);
% fprintf('set: Tprobe = %g us \n',Tprobe/1e-6);
% fprintf('set: M = %g \n',M);
% fprintf('set: Diffusion sensor length = %gx%g um \n',diffusion_sensor_size/mum,diffusion_sensor_size/mum);
% fprintf('Diffusion sensor volume = %g cm3 \n',Vsensor/cm^3);
% fprintf('set: I0 = %g mW/cm2 \n',I0/mW*cm^2);
% fprintf('I_transmitted = %g mW/cm2 \n',I_transmitted/mW*cm^2);
% fprintf('set: optical pumping efficiency = %g \n',optical_pumping);

% counts, shot noise per camera pixel
% Nph_pixel = I_transmitted*Aeff*Tprobe/hbar/om0; %# of photons / pixel / probe pulse
% Nel_pixel = QE*Nph_pixel; % number of electrons / pixel / probe pulse
% delta_Nel_pixel = sqrt(Nel_pixel); % shot noise in the electron number
% Ncounts_pixel = gain*Nel_pixel;

% counts, shot noise per diffusion limited sensor
Nph_sensor = I_transmitted*Adiffusion*Tprobe/hbar/om0; %# of photons / sensor / probe pulse
Nel_sensor = QE*Nph_sensor; % number of electrons / sensor / probe pulse
delta_Nel_sensor = sqrt(Nel_sensor); % shot noise in the electron number
Ncounts = gain*Nel_sensor;
delta_Ncounts = gain*delta_Nel_sensor;
SN_raw = 1/delta_Ncounts;

% fprintf('\n');
% fprintf('PHOTON SHOT NOISE: \n');
% fprintf('Nphotons/pixel = %g \n',Nph_pixel);
% fprintf('Ncounts/pixel = %g \n',Ncounts_pixel);
% fprintf('Ncounts/sensor = %g  \n',Ncounts);
% fprintf('Ncounts/sensor, rel. error = %0.2g \n',SN_raw);
% 
% fprintf('\n');
% fprintf('The rest of the values are for a diffusion-limited sensor \n');

%% OD error
% image and reference typically have counts of the same order of magnitufe
% for simplicity, will assume they are equal. OD should then just be zero
% OD = -log(Nimage/Nref);
% delta_OD = sqrt( (dOD/dNi * deltaNi)^2 + (dOD/dNref * deltaNref)^2 )

Ni = Ncounts;
Nref = Ncounts;
delta_Ni = delta_Ncounts;
delta_Nref = delta_Ncounts;

dODdNi = -1/Ni; %partial derivative OD wrt Ni
dODdNref = 1/Nref; %partial derivative OD wrt Nref

delta_OD = sqrt( (dODdNi * delta_Ni)^2 + (dODdNref * delta_Nref)^2 );

Nshots=150;
Nruns=199;
delta_OD_shot_average=delta_OD/sqrt(Nshots);
delta_OD_run_average=delta_OD/sqrt(Nshots)/sqrt(Nruns);

% fprintf('delta_OD = %0.2g (single shot) \n',delta_OD);
% fprintf('delta_OD = %0.2g (%g averaged shots) \n',delta_OD_shot_average,Nshots);
% fprintf('delta_OD = %0.2g (%g averaged %g-shot runs) \n',delta_OD_run_average,Nruns,Nshots);

%% shot noise limited sensitivity (based off Boehi APL 2012)
% deltaODmin = sqrt(hbar*om0/I0 * A * dt_im)
% Omega_min = dt_im/2 * sqrt(deltaODmin/OD)

dt_mw=T1;
T2 = T1; % coherence time
run_length=30; %time for a single run, in seconds

% standard sequence (no pi/2 preparation)
Omega_min = 2/dt_mw * sqrt(delta_OD/(OD87*optical_pumping/3)) * exp(dt_mw/T2);
Bmw_min = hbar/muB*Omega_min;
% averaged shots
Omega_min_average = 2/dt_mw * sqrt(delta_OD_shot_average/(OD87*optical_pumping/3)) * exp(dt_mw/T2);
Bmw_min_average = hbar/muB*Omega_min_average;
Sensitivity = Bmw_min_average*sqrt(run_length);
% averaged runs
Omega_min_runs = Omega_min_average/sqrt(Nruns);
Bmw_min_runs = Bmw_min_average/sqrt(Nruns);
Sensitivity_runs = Bmw_min_runs*sqrt(run_length*Nruns);

% optimised sequence (preparation with pi/2 pulse)
Omega_min_opt = 2/dt_mw *delta_OD/(OD87*optical_pumping/3) * exp(dt_mw/T2);
Bmw_min_opt = hbar/muB*Omega_min_opt;
% averaged shots
Omega_min_average_opt = 2/dt_mw *delta_OD_shot_average/(OD87*optical_pumping/3) * exp(dt_mw/T2);
Bmw_min_average_opt = hbar/muB*Omega_min_average_opt;
Sensitivity_opt = Bmw_min_average_opt*sqrt(run_length);
% averaged runs
Omega_min_runs_opt = Omega_min_average_opt/sqrt(Nruns);
Bmw_min_runs_opt = Bmw_min_average_opt/sqrt(Nruns);
Sensitivity_runs_opt = Bmw_min_runs_opt*sqrt(run_length*Nruns);

% % sensitivity scales to V=1cm3
% V=50*50*100*mum^3;
% Sensitivity_scaled = Sensitivity*sqrt(V);
% Sensitivity_average_scaled = Sensitivity_average*sqrt(V);

% fprintf('\n');
% fprintf('set: dt_mw = %g us \n',dt_mw/mus);
% fprintf('set: run length = %g s \n',run_length);
% fprintf('\n');
% fprintf('Standard sequence (no pi/2 preparation): \n');
% fprintf('Omega_min = 2pi*%0.3g kHz (single shot), 2pi*%0.3g kHz (%g averaged shots), 2pi*%0.3g kHz (%g averaged %g-shot runs)  \n',...
%     Omega_min/2/pi/kHz,Omega_min_average/2/pi/kHz,Nshots,Omega_min_runs/2/pi/kHz,Nruns,Nshots);
% fprintf('Bmw_min = %0.3g nT (single shot), %0.3g nT (%g averaged shots), %0.3g nT (%g averaged %g-shot runs)\n',...
%     Bmw_min/1e-9,Bmw_min_average/1e-9,Nshots,Bmw_min_runs/1e-9,Nruns,Nshots);
% fprintf('Sensitivity = %0.3g nT/sqrt(Hz) (%g averaged shots), %0.3g nT/sqrt(Hz) (%g averaged %g-shot runs)\n',...
%     Sensitivity/1e-9,Nshots,Sensitivity_runs/1e-9,Nruns,Nshots);
% fprintf('\n');
% fprintf('Optimised sequence (preparation with pi/2 pulse): \n');
% fprintf('Omega_min = 2pi*%0.3g kHz (single shot), 2pi*%0.3g kHz (%g averaged shots), 2pi*%0.3g kHz (%g averaged %g-shot runs)  \n',...
%     Omega_min_opt/2/pi/kHz,Omega_min_average_opt/2/pi/kHz,Nshots,Omega_min_runs_opt/2/pi/kHz,Nruns,Nshots);
% fprintf('Bmw_min = %0.3g nT (single shot), %0.3g nT (%g averaged shots), %0.3g nT (%g averaged %g-shot runs)\n',...
%     Bmw_min_opt/1e-9,Bmw_min_average_opt/1e-9,Nshots,Bmw_min_runs_opt/1e-9,Nruns,Nshots);
% fprintf('Sensitivity = %0.3g nT/sqrt(Hz) (%g averaged shots), %0.3g nT/sqrt(Hz) (%g averaged %g-shot runs)\n',...
%     Sensitivity_opt/1e-9,Nshots,Sensitivity_runs_opt/1e-9,Nruns,Nshots);

%% Atomic Projection Noise

Tmeas = 1; % coherence time
atoms_interacting=1/3; % maximum fraction of prepared atoms that can interact with the mw (1/3 of F=1 atoms)
Nat = density_87Rb*atoms_interacting*Vsensor;
dBatom = hbar/muB/sqrt(Nat*T2*Tmeas);

% fprintf('\n');
% fprintf('ATOMIC PROJECTION NOISE: \n');
% fprintf('\n');
% fprintf('set: T2 = %g mus \n',T2/mus);
% fprintf('set: Tmeas = %g s \n',Tmeas);
% fprintf('Nat = %0.3g atoms per sensor \n',Nat);
% fprintf('dBmw = %0.3g fT/sqrt(Hz) \n',dBatom/fT);

% ratio, optimised shot-noise sensitivity to atomic projection noise:
ratio=Sensitivity/dBatom;
ratio_opt=Sensitivity_opt/dBatom;
% fprintf('\n');
% fprintf('Shot noise is %0.2g times atomic projection noise in the currenct setup, \n',ratio);
% fprintf('or %0.2g times with an optimised sequence \n',ratio_opt);


output.T1=T1;
output.T2=T2;
output.delta_OD=delta_OD;
output.delta_OD_shot_average=delta_OD_shot_average;
output.Omega_min=Omega_min;
output.Omega_min_opt=Omega_min_opt;
output.Sensitivity=Sensitivity;
output.Sensitivity_opt=Sensitivity_opt;
output.Nat=Nat;
output.dBatom=dBatom;
output.diff_coeff=diff_coeff;
output.diffusion_sensor_size=diffusion_sensor_size;
output.delta_OD_run_average=delta_OD_run_average;
output.Ncounts=Ncounts;
output.delta_Ncounts=delta_Ncounts;
output.OD=OD;
output.OD87=OD87;
end