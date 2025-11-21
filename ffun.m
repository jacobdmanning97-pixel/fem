function [localmat] = ffun(xy_pts, triag_no) 
%
% This function computes, the values for ffun 
% at the requested xy_pts points in triangle triag_no.
%  The vector of values is returned in localmat.
%  
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type
global quad_rul num

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%
npts = size(xy_pts,1) ;
x = xy_pts(:,1)';
y = xy_pts(:,2)';

%% localmat is a vector of values
%'num' is used to change between the two problems in the homework
if num == 1
    localmat = 2*pi*x.*y+10;
else
    localmat = 2*(x.*sin(2*pi*x.*y)+x)+3*(sin(2*pi*x.*y)+2*pi*x.*y.*cos(2*pi*x.*y)+1)+4*pi^2*x.*y.^2.*sin(2*pi*x.*y)+ ...
                4*pi^2*x.^3.*sin(2*pi*x.*y)-4*pi*y.*cos(2*pi*x.*y)+4*pi*x.^2.*cos(2*pi*x.*y);
end





