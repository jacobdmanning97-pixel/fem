function [localmat] = Vfun(xy_pts, triag_no) 
%
% This function computes, the values for Vfun -- the
% velocity function
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
xpts = xy_pts(:,1) ;
ypts = xy_pts(:,2) ;

localmat = zeros(2, npts) ;

%% localmat is a vector of values here the function is [-x , y]
if num ==1 
    localmat(1,:) = -xpts.' ;
    localmat(2,:) = ypts.' ;
else
    localmat(1,:) = 3*ones(1,npts) ;
    localmat(2,:) = 2*ones(1,npts)  ;
end
