function [quad_pts, quad_wghts] = quad_73
%
% This function contains weights and quadrature points
% which are exact for polynomials of degree three.
%

quad_pts(1,:) = [0.5, 0.0] ;
quad_pts(2,:) = [0.5, 0.5] ;
quad_pts(3,:) = [0.0, 0.5] ;
quad_pts(4,:) = [0.0, 0.0] ;
quad_pts(5,:) = [1.0, 0.0] ;
quad_pts(6,:) = [0.0, 1.0] ;
quad_pts(7,:) = [1/3, 1/3] ;


quad_wghts(1:3) = 8/120* ones(1,3) ;
quad_wghts(4:6) = 3/120* ones(1,3) ;
quad_wghts(7) = 27/120 ;
