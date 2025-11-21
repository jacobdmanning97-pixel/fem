function [quad_pts, quad_wghts] = quad_32
%
% This function contains weights and quadrature points
% which are exact for polynomials of degree 2.
%

quad_pts(1,:) = [0.5, 0.0] ;
quad_pts(2,:) = [0.5, 0.5] ;
quad_pts(3,:) = [0.0, 0.5] ;

quad_wghts = [1/6 ; 1/6; 1/6] ;

