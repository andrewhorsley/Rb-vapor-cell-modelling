function output = BufferShiftBroadeningFunc(T,Pbuffer_fill,Tfill,bgas)
% Outputs line shifts and broadening for 87Rb optical (D1/2) and microwave
% (hyperfine ground state) transitions. Also outputs the Rb diffusion
% constant and mean-free-path. The script can be adapted for other 
% alkali species by inputting the appropriate constants for that species
% 
% see Andrew Horsley PhD thesis, Ch.2, for detailed theory: 
% High Resolution Field Imaging with Atomic Vapor Cells, Basel 2015
%
% Andrew Horsley 2017

GHz=1e9;
clight=2.99792458e8; % Speed of Light
kB=1.3806504e-23; % Boltzmann's Constant
m87Rb=1.443160648e-25; % Rb87 Atomic Mass
om0=2*pi*384.2304844685e12; % Frequency D2 Transition (also use for D1 line broadening - error is ~2%)
p_atmoshperic=101325; % atmospheric pressure, in Pascals

BufferGasProperties;

n_buffer= Pbuffer_fill/kB/Tfill; % derived buffer density, in m^3
v_av_rb= sqrt(8*kB.*T/pi/m87Rb); %Average Rb velocity

for i=1:length(Pbuffer_fill)
    t_bgas=bgas(i,:);
    
    v_av_bg(i)=sqrt(8*kB*T/pi/mass.(t_bgas));
    vrel_rbbg(i) = sqrt(v_av_rb^2+v_av_bg(i)^2);
    t_gamma_rbbg_total(i)=n_buffer(i)*sigma_total.(t_bgas)(1)*vrel_rbbg(i)*(T/sigma_total.(t_bgas)(2))^sigma_total.(t_bgas)(3); % buffer T1 relaxation rate
    t_gamma_1bg(i)=n_buffer(i)*sigma_1.(t_bgas)(1)*vrel_rbbg(i)*(T/sigma_1.(t_bgas)(2))^sigma_1.(t_bgas)(3); % buffer T1 relaxation rate
    t_gamma_2bg(i)=n_buffer(i)*sigma_2.(t_bgas)(1)*vrel_rbbg(i)*(T/sigma_2.(t_bgas)(2))^sigma_2.(t_bgas)(3); % buffer T2 relaxation rate

    t_RQ(i) = n_buffer(i)*sigmaQ.(t_bgas)(2)*vrel_rbbg(i); % quenching rate for D2 line
    
    t_OptBroad(i)=opt_broad.(t_bgas)(1)*Pbuffer_fill(i)*(T/opt_broad.(t_bgas)(2))^opt_broad.(t_bgas)(3)*(opt_broad.(t_bgas)(2)/Tfill); %optical buffer broadening
    t_OptShift(i)=opt_shift.(t_bgas)(1)*Pbuffer_fill(i)*(T/opt_shift.(t_bgas)(2))^opt_shift.(t_bgas)(3)*(opt_shift.(t_bgas)(2)/Tfill); %optical line shift
    
    t_mw_shift(i)= Pbuffer_fill(i)*(mw_shift.(t_bgas)(4)/Tfill)*...
        (mw_shift.(t_bgas)(1) + ...
        mw_shift.(t_bgas)(2)*(T-mw_shift.(t_bgas)(4)) + ...
        mw_shift.(t_bgas)(3)*(T-mw_shift.(t_bgas)(4))^2); % hyperfine line shift
    
    t_D(i) =D0.(t_bgas)(1)*p_atmoshperic/ (Pbuffer_fill(i)*(D0.(t_bgas)(2)/Tfill)) *(T/D0.(t_bgas)(2))^D0.(t_bgas)(3); %Diffusion constant D, in m^2/s
    
    t_mfp_D0(i) = 3*t_D(i)/vrel_rbbg(i); % mean free path, from D0
    output.t_mfp_sigma(i) = 1/(n_buffer(i)*sigma_total.(t_bgas)(1)*sqrt(1+mass.(t_bgas)/m87Rb)); % mean free path, from sigma_rbbg_total   

end

% combine properties of the (one or many) buffer gases as single quantities
output.gamma_rbbg_total=sum(t_gamma_rbbg_total); % total Rb-buffer collision rate
output.RQ = sum(t_RQ); % quenching rate
output.OptBroad=sum(t_OptBroad); 
output.OptShift=sum(t_OptShift);
output.MW_shift=sum(t_mw_shift);
output.mfp_rbbg =1/sum(1./output.t_mfp_sigma);

output.diff_coeff = (sum(1./t_D))^-1;

% doppler linewidth
output.doppler_opt=om0/2/pi/clight *sqrt(8*kB*T*log(2)/m87Rb);
output.doppler_mw=6.835*GHz/clight *sqrt(8*kB*T*log(2)/m87Rb);
% dicke narrowed linewidths
output.dicke_opt=4*pi*output.diff_coeff/(clight/(om0/2/pi))^2;
output.dicke_mw=4*pi*output.diff_coeff/(clight/(6.835*GHz))^2;


