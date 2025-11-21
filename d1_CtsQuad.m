function [CQuadVal, DerivCQuadVal] = d1_CtsQuad(quad_pts)
%
% This function computes the values of the continuous
% quadratic basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference interval (0, 1).
%

CQuadVal(1,:) = 2.0*(quad_pts-1/2).*(quad_pts-1);
CQuadVal(2,:) = -4.0*quad_pts.*(quad_pts-1);
CQuadVal(3,:) = 2.0*quad_pts.*(quad_pts-1/2);


DerivCQuadVal(1,:) = 2.0*(quad_pts-1/2)+2.0*(quad_pts-1) ;
DerivCQuadVal(2,:) = -4.0*quad_pts-4.0*(quad_pts-1)  ;
DerivCQuadVal(3,:) = 2.0*quad_pts+2.0*(quad_pts-1/2)  ;
