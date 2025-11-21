function [localmat] = gfun(xy_pts) 
%
% This function computes, the values for Gfun -- the
% Dirichlet boundary data
% at the requested xy_pts points 
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
x_pts = xy_pts(:,1) ;
y_pts = xy_pts(:,2) ;

%% localmat is a vector of values
if num == 1
    localmat = 2* x_pts.^2 + pi* x_pts.*y_pts + 7 ;
else
    localmat = (x_pts .*sin(2*pi *x_pts .*y_pts) + x_pts ).' ;
end




