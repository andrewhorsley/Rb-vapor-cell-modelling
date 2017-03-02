% Computes various properties for a vapor cell, and prints these in the
% command window
%
% Andrew Horsley 2017

clear all
clc
constants; units; 

%% Input parameters

%%%%% Cell temperature %%%%%%%
T=60 +273.15; %cell temperature, in Kelvin

%%%%% Rb properties %%%%%%%
rb87_abundance = 0.98; %natural abundance of 87Rb. Set to 1 for 100% 87Rb cell
Rb=RbProperties(T,rb87_abundance);

%%%%% Cell and laser beam dimensions
cell_shape = 'cylinder'; % 'sphere', 'cylinder', rectangle'
if strcmp(cell_shape, 'sphere')
    cell_radius = 0.5*mm; % (optical path length)
    cell_dimensions = cell_radius;
    optical_path_length = cell_radius;
elseif strcmp(cell_shape, 'cylinder')
    cell_radius = 0.5*mm; 
    cell_thickness = 37.5*mm; % (optical path length
    cell_dimensions = [cell_radius, cell_thickness];
    optical_path_length = cell_thickness;
elseif strcmp(cell_shape, 'rectangle')
    a = 6*mm; % transverse dimensions
    b = 6*mm; % transverse dimensions
    c = 200*mum; % cell thickness (optical path length)
    cell_dimensions = [a, b, c];
    optical_path_length=c;
end

%%%%% Buffer gas properties %%%%%%%
bgas=['N2']; % buffer gas types (as many as you want, e.g. ['Kr'; 'N2'])
Pbuffer_fill=[11]*torr; %buffer gas pressure, e.g. [10; 25]*mbar. Vector length must match bgas
Tfill = 25+273.15; % temperature at which Pbuffer is defined, default 25+273.15 kelvin

%%%% alternative input method for buffer gas pressure:
% P0total = 105 *mbar; % for matching pbuffer to observed T1
% P1_frac = 0.8; P2_frac = 1-P1_frac;
% P0_1 = P1_frac*P0total; P0_2 = P2_frac*P0total;
% Pbuffer_fill=[P0_1 P0_2];
% Tfill = 25+273.15; % temperature at which Pbuffer is defined, default 25+273.15 kelvin

%%%%% Laser frequency and detuning (Delta) %%%%%%%
Dline = 'D2'; % options are 'D1' or 'D2'

% Delta=1.28e9; %for max OD with zero buffer - centre of doppler-broadened D2 line
Delta = 2*pi*(266.65/2)*1e6; %laser locked to Fg=2->Fe=2/3 crossover
% Delta = Delta - 2*pi*80e6 + 2*pi*120e6; % to apply a detuning offset. otherwise comment out
 

%% Collisional broadening and shifts

Buffer = BufferShiftBroadeningFunc(T,Pbuffer_fill,Tfill,bgas,Rb);

%% T1 and T2 lifetimes

[T1, T2, gamma] = ...
    RelaxationFunc(cell_shape,cell_dimensions,T,Pbuffer_fill,Tfill,bgas,Rb);
% gamma.SE1, gamma.SE2, gamma.walls, gamma.bg1, gamma.bg2

%% Diffusion rates:
diff_time.a = 10*mus;
diff_time.b = 20*mus;
diff_time.c = 50*mus;

diff_dist.a = DiffusionTimeFunc(diff_time.a,Buffer.diff_coeff);
diff_dist.b = DiffusionTimeFunc(diff_time.b,Buffer.diff_coeff);
diff_dist.c = DiffusionTimeFunc(diff_time.c,Buffer.diff_coeff);
diff_dist.T1 = DiffusionTimeFunc(T1,Buffer.diff_coeff);

%% Optical depth:
n_buffer= Pbuffer_fill/kB/Tfill; % derived buffer density, in m^3

Esus_out=Esus_func(Dline,0,T,optical_path_length,rb87_abundance,Delta,n_buffer,bgas);
OD=Esus_out.OD_p;
OD85=Esus_out.OD85_p;
OD87=Esus_out.OD87_p;

%% Printing out results:

fprintf('/MATLAB/Vapor Cell (Matlab/Hellma_cell/Hellma_100_cell_data.m \n');
fprintf('\n');
fprintf('Cell Temperature = %g degC \n', T-273.15);
fprintf('Optical path length = %g mm \n', optical_path_length/mm);
for i=1:length(Pbuffer_fill)
fprintf('%s %g = %0.3g mbar %s (fill)  \n', 'Buffer gas',i, Pbuffer_fill(i)/mbar, bgas(i,:));
end
fprintf('MW shift = %0.3g kHz \n', Buffer.MW_shift/kHz);
fprintf('Optical buffer shift = %0.3g MHz \n', Buffer.OptShift/MHz);
fprintf('Optical buffer broadening = %0.3g GHz \n', Buffer.OptBroad/GHz);
fprintf('Doppler optical = %0.3g GHz \n', Buffer.doppler_opt/GHz);
fprintf('Doppler mw = %0.3g kHz, Dicke narrowed = %0.3g Hz \n', Buffer.doppler_mw/kHz,Buffer.dicke_mw);
fprintf('OD = %0.4g (OD87 = %0.3g, OD85 = %0.3g) \n', OD, OD87, OD85);
fprintf('Rb87 vap partial pressure = %0.3g mbar \n', Rb.P/mbar*rb87_abundance);
fprintf('Total Rb density = %0.3g atoms/cm^3 \n', Rb.density_total*cm^3);
fprintf('87Rb density = %0.3g atoms/cm^3 (87Rb abundance = %g) \n', Rb.density_87*cm^3, rb87_abundance);
fprintf('\n');
fprintf('Mean free path, Rb-Rb = %0.3g mm \n', Rb.lambda_rbrb/mm);
for i=1:length(Pbuffer_fill)
fprintf('Mean free path, Rb-buffer %g = %0.3g nm\n',i, Buffer.t_mfp_sigma(i)/nm);
end
fprintf('Mean free path, Rb-buffer total = %0.3g um\n', Buffer.mfp_rbbg);
fprintf('Collision rate Rb-Rb (total) = %g /s \n', Rb.gamma_rbrb_total);
fprintf('Collision rate Rb-buffer (total) = %g /s \n', Buffer.gamma_rbbg_total);
fprintf('quenching rate = %g /s \n', Buffer.RQ);
fprintf('Diffusion constant = %g cm^2/s \n', Buffer.diff_coeff/cm^2);
fprintf('Diff. dist. (um) %gus=%0.3g, %gus=%0.3g, %gus=%0.3g \n',... 
    diff_time.a/mus,diff_dist.a/mum, ...
    diff_time.b/mus, diff_dist.b/mum, ...
    diff_time.c/mus, diff_dist.c/mum);
fprintf('Diff. dist. (um) during T1 = %0.3g \n', diff_dist.T1/mum);
fprintf('\n');
fprintf('T1 = %g us \n', T1/mus); %Based off 1/T1 = 1/TG+1/Tex, Arditi 1964
fprintf('T2 = %g us \n', T2/mus); 
fprintf('gamma.SE1 = %0.3g /s \n', gamma.SE1);
fprintf('gamma.SE2 = %0.3g /s \n', gamma.SE2);
fprintf('gamma.walls = %0.3g /s \n', gamma.walls);
fprintf('gamma.bg1 = %g /s \n', gamma.bg1);
fprintf('gamma.bg2 = %g /s \n', gamma.bg2);


