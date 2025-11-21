function [localmat] = efun(x_pts, isub) 
%
% This function represents the e function in the 
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
localmat = 5*x_pts + 2  ;


