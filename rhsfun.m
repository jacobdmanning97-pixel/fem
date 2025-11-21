function [localmat] = rhsfun(x_pts, isub, num) 
%
% This function represents the rhs function in the 
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
    localmat = 3*x_pts.^2+(3*x_pts-3).*(5*x_pts+2)-6;
elseif num == 2 
    %3x^2-3
    localmat = 6*x_pts.^3+(3*x_pts.^2-3).*(5*x_pts+2)-6*(2*x_pts+1)-12*x_pts;
else
    %x^4+1
    localmat = 4*x_pts.^5+(x_pts.^4+1).*(5*x_pts+2)-8*x_pts.^3-12*x_pts.^2.*(2*x_pts+1)  ;
end
