%%%%%%%%%%%%%%%%%%%
%
% This files contains the command for computing the L2 error
% for u and Gradu and the H1 error for the u
%
% Note: If utruefun, Gradutruefun and are set to zero this
% routine computes the appropriate norms for u.

function [uL2Error, GraduL2Error, H1uError]=CalcErr(nodeco, elnode)

quad_pts = [ ] ;
quad_wghts = [ ] ;


GraduErr = zeros(2,1) ;
uErr = 0.0 ;

for itrg = 1:size(elnode,1)
      
   % Description of triangle.
	cotri(1:3,1) = nodeco(elnode(itrg, 1:3), 1) ;
	cotri(1:3,2) = nodeco(elnode(itrg, 1:3), 2) ;
    
	Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
        	(cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
	detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
	JInv = inv(Jmat) ;


	% Evaluation of quadrature points and quadrature weights.
%	[quad_pts, quad_wghts] = feval('quad_73') ;
	[quad_pts, quad_wghts] = feval('quad_137') ;
	nqpts = size(quad_pts,1) ;

	% Adjust points and weights to account for size of true triangle.
	xy_pts = ( Jmat * quad_pts.' ).' ;
	xy_pts(:,1) = cotri(1,1) + xy_pts(:,1) ;
	xy_pts(:,2) = cotri(1,2) + xy_pts(:,2) ;
	quad_wghts = detJ * quad_wghts ;

	% Evaluate U, and grad U  at the quadrature points. 
        Gradu_vals = GradUFun(xy_pts, itrg) ;
        ufun_vals = UFun(xy_pts, itrg) ;
   
        Gradutru_vals = Gradutrue(xy_pts, itrg) ;
        utru_vals = utruefun(xy_pts, itrg) ;
   
        uErr = uErr + quad_wghts * ( (utru_vals - ufun_vals).^2 ).' ;

        temp(1:nqpts) = (Gradutru_vals(1,:) - Gradu_vals(1,:) ).^2 ;
        GraduErr(1,1) = GraduErr(1,1) +  quad_wghts * temp.'  ;
        temp(1:nqpts) = (Gradutru_vals(2,:) - Gradu_vals(2,:) ).^2 ;
        GraduErr(2,1) = GraduErr(2,1) + quad_wghts * temp.'  ;


end

%%%

uL2Error = sqrt(uErr);
   
GraduL2Error = sqrt( GraduErr(1,1) + GraduErr(2,1) );

H1uError = sqrt( uErr + GraduErr(1,1) + GraduErr(2,1) );
