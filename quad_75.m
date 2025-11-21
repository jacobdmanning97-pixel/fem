function [quad_pts, quad_wghts] = quad_75
%
% This function contains weights and quadrature points
% which are exact for polynomials of degree five.
%

quad_pts(1,:) = [(6+sqrt(15))/21, (9-2*sqrt(15))/21] ;
quad_pts(2,:) = [(6+sqrt(15))/21, (6+sqrt(15))/21] ;
quad_pts(3,:) = [(9-2*sqrt(15))/21, (6+sqrt(15))/21] ;
quad_pts(4,:) = [(6-sqrt(15))/21, (6-sqrt(15))/21] ;
quad_pts(5,:) = [(9+2*sqrt(15))/21, (6-sqrt(15))/21] ;
quad_pts(6,:) = [(6-sqrt(15))/21, (9+2*sqrt(15))/21] ;
quad_pts(7,:) = [1/3, 1/3] ;


quad_wghts(1:3) = (155+sqrt(15))/2400* ones(1,3) ;
quad_wghts(4:6) = (155-sqrt(15))/2400* ones(1,3) ;
quad_wghts(7) = 270/2400 ;
