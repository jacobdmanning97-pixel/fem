function [localmat] = uTrue(x_pts, isub, num) 
%
% This function represents the diffusive function in the 
% differential equation.
%  
%  
%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global xpts nnds
global Global_r  Global_s  Global_u
global rad_bas_type  str_bas_type  vel_bas_type
global quad_rul
%
%

nevalpts = size(x_pts,1) ;

if num == 1
    %3x-3
    localmat = 3*x_pts-3 ;
elseif num == 2 
    %3x^2-3
    localmat = 3*x_pts.^2-3;
else
    %x^4+1
    localmat = x_pts.^4+1 ;
end
