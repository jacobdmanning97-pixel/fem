function [val] = d1_ip_fun(ab, quadrule, func)
%
% This function approximates the integral of the function func
% over the interval (ab(1), ab(2)) using the quadrature rule quadrule.
%  

% Description of subinterval.
xleft = ab(1)  ;
xright = ab(2)  ;
hsub = xright - xleft ;

% Evaluation of quadrature points and quadrature weights.
[quad_pts, quad_wghts] = feval(quadrule) ;
nqpts = size(quad_pts,1) ;

% Adjust points and weights to account for size of the interval.
quad_pts = xleft + hsub* quad_pts ;
quad_wghts = hsub * quad_wghts ;

% Evaluate the function at the quadrature points.
fun_vals = feval(func, quad_pts) ;

% Now to approximate the integral
val = fun_vals * quad_wghts ;


