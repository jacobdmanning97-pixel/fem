function [localmat] = utruefun(xy_pts, triag_no)  
%
% This function represents the true solution
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
xpts = xy_pts(:,1) ;
ypts = xy_pts(:,2) ;

%% localmat is a vector of values
if num == 1
    localmat = (2* xpts.^2 + pi* xpts.*ypts + 7).' ;
else
    localmat = (xpts .*sin(2*pi *xpts .*ypts) + xpts ).' ;
end



