%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%
% This is a driver routine for the Finite Element approximation of:
%   -Grad . Kfun Grad u + Vfun Grad u + qfun u = ffun  , for x in R,
%      subject to   u = gfun  on the boundary   of R
%
%  Kfun may either denote a scalar or a 2x2 matrix.

function [uL2Error, GraduL2Error, H1uError] = Drv_ConDiff()

%% In the implemenation the unknown u is associated with GlobalV

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type
global quad_rul num mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Define the problem and the solution method
%%

% We begin by specifying the type of approximation elements
% we are using for the u

% vel_bas_type = 'CtsQuad' ;
% quad_rul = 'quad_75';
% mesh = 'R12x12';
% num = 1;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Set up the Geometry and the Boundary/Initial conditions.

% We assume that the subroutine will pass back 
% nodeco, elnode, and bdynde.
% nodeco will contain the x and y coordinate of the nodes
% elnode will contain the definition of the triangulation
% bdynde will contain the definition of the boundary in a
% counter clockwise direction. 

% Enter the grid file
%%
% % % [nodeco, elnode, eldegr, elnbgh, bdynde, ellevels, curels] =  ...
% % %           ADgridrec([-1 ; 1], [2], [0 ; 1], [2]) ;
% % % 
% % % % The routine midEDGEgen is then called to generate the midedges
% % % % for the mesh
% % % [nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;

% The mesh file already had the midedge etc. created in the appropiate format.
if strcmp(mesh, 'R4x4') == 1
    R4x4 ;
elseif strcmp(mesh, 'R6x6') == 1
    R6x6;
elseif strcmp(mesh, 'R8x8') == 1
    R8x8;
elseif strcmp(mesh, 'R10x10') == 1
    R10x10;
elseif strcmp(mesh, 'R12x12') == 1
    R12x12;
end

%% At this stage you should have the arrays:
%% nodeco, elnode, bdynde, bdyedge, nVert, nedge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Set up and initialize the Global Solution vector.
%

if strcmp(vel_bas_type, 'CtsQuad') == 1 
   GlobalV = zeros((nVert + nedge), 1) ;
        
elseif  strcmp(vel_bas_type, 'CtsLin') == 1 
   GlobalV = zeros(nVert, 2) ;
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4. Set up the approximating linear system.
%
% Set up the approximating linear system.

%% Firstly we allocate the space for the linear system and set up the local to
%% Guess an allocation for the coefficient matrix and RHS
Sdim = size(GlobalV,1) ;
Acoeff = spalloc(Sdim, Sdim, 9*Sdim) ;
RHSvec = zeros(Sdim,1) ;

ntri = size(elnode,1) ;

for itrg = 1:ntri
      
      Alocal = [ ] ;
      RHSloc = [ ] ;

      % Set up unknown solution vector mapping to global velocity vector
      if strcmp(vel_bas_type, 'CtsQuad') == 1
            Vstart = [elnode(itrg,1:3) , nVert + elnode(itrg,4:6) ];
            GlTrgVe =  Vstart.' ;
            
      elseif strcmp(vel_bas_type, 'CtsLin') == 1
            Vstart = elnode(itrg,1:3);
            GlTrgVe = Vstart.' ;
      end

      

      %%%%%%%%%%%%%%%%
      % Evaluate the integrals.

       Alocal = inner_prod_Grad_ten0_Grad_ten0(itrg, quad_rul, ...
          'kfun', vel_bas_type, vel_bas_type) ;
       Alocal = Alocal + inner_prod_ten0_Grad_ten0_Vec(itrg, quad_rul, ...
          'Vfun', vel_bas_type, vel_bas_type) ;
       Alocal = Alocal  + inner_prod_ten0_ten0(itrg, quad_rul, ...
          'qfun', vel_bas_type, vel_bas_type) ;

      RHSloc = inner_prod_ten0(itrg, quad_rul, ...
          'ffun', vel_bas_type) ;
         
      %% Distribute the matrix and the RHS vector to the system. 
      Acoeff(GlTrgVe , GlTrgVe) = Acoeff(GlTrgVe , GlTrgVe) + Alocal ;
      RHSvec(GlTrgVe) = RHSvec(GlTrgVe) + RHSloc ;
   
end

%%%%%%%%%%%%%%%%
% Now apply the boundary conditions.
%%% I am assuming here that Dirichlet BC apply 
    nbdy = size(bdynde,1) ;
    
    for ibdy = 1:nbdy
          eqn = bdynde(ibdy,1) ;
          val = gfun([nodeco(eqn,1), nodeco(eqn,2)]) ;
          RHSvec = RHSvec - val*Acoeff(:,eqn) ;
          Acoeff(:, eqn) = 0 ;
          Acoeff(eqn, :) = 0 ;
          Acoeff(eqn, eqn) = 1.0 ;
          RHSvec(eqn) = val ;
    end

%%% For quadratics also have to take care of the unknown along the boundary 
%%% edges. Left for students!

if strcmp(vel_bas_type, 'CtsQuad') == 1

   nbdye = size(bdyedge,1) ;

   for ibdye = 1:nbdye
      eqn = nVert + bdyedge(ibdye,1) ;
      val = gfun([nodeco(eqn,1), nodeco(eqn,2)]) ;
      RHSvec = RHSvec - val*Acoeff(:,eqn) ;
      Acoeff(:, eqn) = 0 ;
      Acoeff(eqn, :) = 0 ;
      Acoeff(eqn, eqn) = 1.0 ;
      RHSvec(eqn) = val ;
   end
end
%%%%
%%%%

%%%%% Scale the linear system and then solve
%% Scale the matrix
MaxVal = max(abs(Acoeff).'); % This give us the max value in each row
scalvec = (1./MaxVal).'; % This gives us our scalar vector
Acoeff = spdiags(scalvec,0,Sdim,Sdim)*Acoeff; 
RHSvec = scalvec.* RHSvec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Now to solve the problem

GlobalV = Acoeff \ RHSvec ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: Postprocessing


%%  Plot the solution 
    %graphCD2d ;

%% Calculate the error
     [uL2Error, GraduL2Error, H1uError] = CalcErr(nodeco, elnode);






      
