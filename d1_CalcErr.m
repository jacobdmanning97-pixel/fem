%%%%%%%%%%%%%%%%%%%
%
% This files contains the command for computing the L2 error
% for u and Du.
%
% Note: If uTrue and DuTrue are set to zero this
% routine computes the appropriate norms for the radial function.

function [uErr, DuErr] = d1_CalcErr(uTrue_Flag,nnds, num)

global quad_rul
global xpts

uErr = 0 ;
DuErr = 0 ;
for isub = 1 : nnds-1
      
   % Description of triangle.
    xleft = xpts(isub) ;
    xright = xpts(isub + 1) ;
    hsub = xright - xleft ;


	% Evaluation of quadrature points and quadrature weights.
	[quad_pts, quad_wghts] = feval(quad_rul) ;
	nqpts = size(quad_pts,1) ;

	% Adjust points and weights to account for size of true subinterval.
    quad_pts = xleft + hsub* quad_pts ;
    quad_wghts = hsub * quad_wghts ;

   % Evaluate the u and its derivative at the quadrature points. 
   u_vals = uFun(quad_pts, isub)' ;
   Du_vals = DuFun(quad_pts, isub)' ;

   
   if uTrue_Flag == true
       u_Truevals = uTrue(quad_pts, isub, num) ;
       Du_Truevals = DuTrue(quad_pts, isub, num) ;
       
       uErr = uErr + ( (u_vals - u_Truevals).^2 )' * quad_wghts  ;
       DuErr = DuErr +  ( (Du_vals - Du_Truevals).^2 )' * quad_wghts ;

   else
       uErr = uErr + ( (u_vals).^2 )' * quad_wghts ;
       DuErr = DuErr +  ( (Du_vals).^2 )' * quad_wghts ; 

   end
end
uErr = sqrt( uErr ) ;
DuErr = sqrt( DuErr ) ;
if uTrue_Flag ==true   
   %disp(['L2-error for u in the approximation ', num2str(uErr)]) ;
   %disp(['L2-error for Du in the approximation ', num2str(DuErr)]) ;
else   
   %disp(['L2-size for u in the approximation ', num2str(uErr)]) ;
   %disp(['L2-size for Du in the approximation ', num2str(DuErr)]) ;
end
