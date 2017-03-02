function diff_dist = diffusion_time_func(time,diffusion_coeff)
% Calculates the diffusion distance during a given time for a given 
% diffusion constant, assuming 2D diffusion. Units are SI
%
% Andrew Horsley 2017

diff_dist = sqrt(2*time*diffusion_coeff);
