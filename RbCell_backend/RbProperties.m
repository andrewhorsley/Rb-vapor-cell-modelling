function output = RbProperties(T,rb87_abundance)
% outputs common properties for Rb vapor at a given temperature and 
% isotopic abundance
%
% P:                Rb vapor (partial) pressure
% density_total:    Rb vapor density (85Rb+87Rb)
% density_85:       85Rb vapor density
% density_87:       87Rb vapor density
% lambda_rbrb:      mean free path for Rb-Rb collisions
% v_av_rb:          Average Rb velocity
% vrel_rbrb:        Average relative velocity of two Rb atoms
% gamma_rbrb_total: total Rb-Rb collision rate
% gamma_rbrb_SE:    rate of Rb-Rb collisions resulting in spin-exchange relaxation
% Gamma_dipole:     Rb D2-line dipole-dipole broadening
%
% Andrew Horsley, 2015

units; constants;

I=3/2; %nuclear spin
sigma_rbrb_total = 1397*10^-16*cm^2; %total cross section for rb-rb collision, Croucher(1969)
sigma_rbrb_SE = 1.9*10^-14*cm^2; %spin exchange cross section for rb, Walter (2002)


if T-273.15<39.3
    P= 10.^(2.881+4.857-4215./T)*torr; % Steck. Total Rb pressure with Rb in solid phase, in Pa
%     P = 133.323 * 10.^(-94.04826-1961.258./T-0.03771687*T+42.57526*log10(T));
else
    P = 10.^(2.881+4.312-4040./T)*torr; % Steck. Total Rb pressure with Rb in liquid phase, in Torr
%     P = 133.323 * 10.^(15.88253-4529.635./T+0.00058663*T-2.99138*log10(T));
end

density_total=P/kB./T; %Total Rb density at T, in m^3
density_85=density_total*(1-rb87_abundance); %85Rb density at T, in m^3
density_87=density_total*rb87_abundance; %87Rb density at T, in m^3

lambda_rbrb = kB.*T/sqrt(2)/sigma_rbrb_total./P; %mean free paths
v_av_rb= sqrt(8*kB.*T/pi/m87Rb); %Average Rb velocity
vrel_rbrb = sqrt(2)*v_av_rb; %Rb-Rb rel velocity
gamma_rbrb_total=density_total*sigma_rbrb_total.*vrel_rbrb; %total collision rate
gamma_rbrb_SE=density_total*sigma_rbrb_SE.*vrel_rbrb; %SE collision rate

Beta2 = 2*pi*1.10*10^-7*cm^3;
Gamma_dipole = Beta2*density_total;  % Rb D2 dipole-dipole broadening

output.P=P;
output.density_total=density_total;
output.density_85=density_85;
output.density_87=density_87;
output.lambda_rbrb=lambda_rbrb;
output.v_av_rb=v_av_rb;
output.vrel_rbrb=vrel_rbrb;
output.gamma_rbrb_total=gamma_rbrb_total;
output.gamma_rbrb_SE=gamma_rbrb_SE;
output.Gamma_dipole=Gamma_dipole;
end



