function [localmat] = qfun(xy_pts, triag_no) 
%
% This function computes, the values for Qfun -- the
% sink/source coefficient
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

%% localmat is a vector of values
localmat = 2*ones(1,npts) ;






