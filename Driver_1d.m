%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%%  In this file we implement the FEM approximation of 
%%  the convection-diffusion problem:
%%  -kfun(x) u'' + bfun(x) u' + efun(x) u = rhsfun(x), alpha < x < beta
%%   u(alpha) = u0,   u(beta) = u1

%%  The functions kfun, bfun, efun are defined in files
%%  
%

function [uErr, DuErr] = Driver_1d(nnds,alpha,beta, uTrue_Flag, num, type)


%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global xpts 
global Global_r  Global_s  Global_u
global rad_bas_type  str_bas_type  vel_bas_type
global quad_rul

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Define the problem 
%%
%% Variables used:  Global_u -- represents the concentration
%%                 

u0 = uTrue(alpha,0,num);
u1 = uTrue(beta,0,num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Define the Approximation parameters

% We assume that the subroutine will pass back xpts a 
% vector of dimension nnds x 1 containing the nodal values.
%

vel_bas_type = type ;        %% used for the concentration

quad_rul = 'd1_quad_5';       %% Decide which numerical quadrature to be used


%% nnds -- represents the number of nodes (Nx := number of subintervals = nnds - 1)

xpts = linspace(alpha, beta, nnds).' ;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Set up and initialize the Global Solution vectors.
%         Also determine the approx bandwidth of the coefficient matrix

if strcmp(vel_bas_type, 'd1_CtsLin') == 1 
   Global_u = zeros(nnds, 1) ;
   xqpts = xpts ;
   bandw = 3 ;
        
elseif  strcmp(vel_bas_type, 'd1_CtsQuad') == 1 
   Global_u = zeros(2*nnds - 1, 1) ;
   xqpts = linspace(alpha, beta, 2*nnds - 1)' ;
   bandw = 5 ;
     
elseif  strcmp(vel_bas_type, 'd1_CtsCub') == 1 
   Global_u = zeros(3*nnds - 2, 1) ;
   xqpts = linspace(alpha, beta, 3*nnds - 2)' ;
   bandw = 7 ;

end


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Assemble the coefficient matrix and rhs vector 

%% Guess an allocation for the coefficient matrix
udim = size(Global_u,1) ;
Amat = spalloc(udim, udim, bandw*udim) ;
bvec = zeros(udim, 1) ;

for isub = 1 : nnds-1 % loop over the subintervals
      
      Alocal = [ ] ;
      RHSloc = [ ] ;

      % Identify the global unknown coefficients
      if strcmp(vel_bas_type, 'd1_CtsLin') == 1
         GlTrg_u = [isub ; isub + 1] ; 
         
      elseif strcmp(vel_bas_type, 'd1_CtsQuad') == 1
         GlTrg_u = [2*isub - 1; 2*isub; 2*isub + 1] ;
          
      elseif strcmp(vel_bas_type, 'd1_CtsCub') == 1
         GlTrg_u = [3*isub - 2; 3*isub - 1; 3*isub; 3*isub + 1] ;
           
      end
      
      %%%%%%%%%%%%%%%%
      % Evaluate the integrals.

      %% The diffusion term
      Alocal =  d1_ip_Dten0_Dten0(isub, 'kfun', vel_bas_type, vel_bas_type)  ;

      %% The transport term
      Alocal =  Alocal + ...
                d1_ip_ten0_Dten0(isub, 'bfun', vel_bas_type, vel_bas_type)  ;
       
      %% The adsorption term
      Alocal =  Alocal + ...
                d1_ip_ten0_ten0(isub, 'efun', vel_bas_type, vel_bas_type)  ;


      %% The rhs term
      RHSloc = d1_ip_ten0(isub, 'rhsfun', vel_bas_type, num)  ;
  
             
      %% Distribute the matrix to the system.
         
      Amat(GlTrg_u , GlTrg_u) = Amat(GlTrg_u , GlTrg_u) + Alocal ;
      bvec(GlTrg_u) = bvec(GlTrg_u) + RHSloc ;
      
end
%full(Amat)
%% Impose the boundary conditions
%% Doing the following preserves the symmetry (if it exists) for the coeff. matrix
bvec = bvec - u0* Amat(: , 1)  - u1* Amat(: , udim);
Amat(: , 1) = zeros(size(Amat(: , 1))) ;
Amat(1 , :) = zeros(size(Amat(1 , :))) ;   
Amat(: , udim) = zeros(size(Amat(: , udim))) ;
Amat(udim , :) = zeros(size(Amat(udim , :))) ;   

Amat(1,1) = 1.0 ;
bvec(1) = u0 ;
Amat(udim, udim) = 1.0 ;
bvec(udim) = u1 ;


%% Scale the matrix
MaxVal = max(abs(Amat).'); % This give us the max value in each row
scalvec = (1./MaxVal).'; % This gives us our scalar vector
Amat = spdiags(scalvec,0,udim,udim)*Amat; 
bvec = scalvec.* bvec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Now to solve the problem

Global_u = Amat \ bvec ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: Postprocessing


%%  Plot the solution 
    %plot(xqpts, Global_u, '-') ;

%% Calculate the error
[uErr, DuErr] = d1_CalcErr(uTrue_Flag,nnds, num);

