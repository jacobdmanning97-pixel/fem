function [CLinVal, DerivCLinVal] = d1_CtsLin(quad_pts)
%
% This function computes the values of the continuous
% linear basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference interval (0, 1).
%

nqpt = size(quad_pts,1) ;

CLinVal(1,:) = 1.0 - quad_pts(1:nqpt)';
CLinVal(2,:) = quad_pts(1:nqpt)' ;


DerivCLinVal(1,:) = -1*ones(1,nqpt) ;
DerivCLinVal(2,:) =  1*ones(1,nqpt) ;
