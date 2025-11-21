function [localmat] = inner_prod_ten0_Grad_ten0_Vec(triag_no, quad_rul, ...
   vc_fun, ten0a_type, ten0b_type) 
%
% This function computes, for triangle triag_no, the integrals
% of V.\del u*v
%


%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type
global quad_rul num

% Description of triangle.
cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
    
Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
        (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
JInv = inv(Jmat) ;

% Evaluation of quadrature points and quadrature weights.
[quad_pts, quad_wghts] = feval(quad_rul) ;
nqpts = size(quad_pts,1) ;

% Adjust points and weights to account for size of true triangle.
xy_pts = ( Jmat * quad_pts.' ).' ;
xy_pts(:,1) = cotri(1,1) + xy_pts(:,1) ;
xy_pts(:,2) = cotri(1,2) + xy_pts(:,2) ;
quad_wghts = detJ * quad_wghts ;

% Evaluate the scalar multiplier at the quadrature points.
vcfun_vals = feval(vc_fun, xy_pts, triag_no) ;

% Evaluate Basis Functions and their Gradients at quad. points.
[ten0a, Gradten0a] = feval(ten0a_type, quad_pts) ;
nbas0a = size(ten0a,1) ;

[ten0b, Gradten0b] = feval(ten0b_type, quad_pts) ;
nbas0b = size(ten0b,1) ;

% Do appropriate multiplies to get the true Gradients.
for iq = 1:nqpts
   Gradtrue1(:,:,iq) = Gradten0b(:,:,iq) * JInv ;
end

% MATLAB will not take the transpose nor do the multiplication of my
% Gradtrue matrices -- Hence we introduce tempM1 and tempM2
tempM1(:,:) = Gradtrue1(:,1,:) ;
tempM2(:,:) = Gradtrue1(:,2,:) ;

% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   tempM1(:,iq) = quad_wghts(iq) * vcfun_vals(1,iq) * tempM1(:,iq) ; 
   tempM2(:,iq) = quad_wghts(iq) * vcfun_vals(2,iq) * tempM2(:,iq) ; 
end

mat1 = ten0a*(tempM1+tempM2)';

localmat = [ mat1 ] ; 
