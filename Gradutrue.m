function [localmat] = Gradutrue(xy_pts, triag_no)  
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

localmat = zeros(2,1,npts) ;

%% localmat is a vector of values
if num == 1
    localmat(1,:) = (4* xpts + pi *ypts).' ;
    localmat(2,:) = ( pi *xpts ).' ;
else
    localmat(1,:) = ( 1 + sin(2*pi* xpts.* ypts) + 2*pi* xpts .*ypts .*cos(2*pi* xpts .*ypts) ).' ;
    localmat(2,:) = ( 2*pi* xpts.^2 .*cos(2*pi* xpts .*ypts) ).' ;
end



