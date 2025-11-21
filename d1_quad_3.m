function [quad_pts, quad_wghts] = d1_quad_3
%
% This function contains weights and quadrature points
% which are exact for polynomials of degree three for evaluting
% the line integral from (0,0) to (1,0).
%

quad_pts = zeros(2,1) ;
quad_wghts = zeros(2,1) ;

quad_pts(1) = (1 - 1/sqrt(3))/2 ;
quad_pts(2) = (1 + 1/sqrt(3))/2 ;


quad_wghts(1:2) = 0.5* ones(1,2) ;

