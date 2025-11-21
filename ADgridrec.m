function [nodeco, elnode, eldegr, elnbgh, bdynde, ellevels, curels] =  ...
         ADgridrec(xvec, nxvec, yvec, nyvec)
%
% This function generates a permutation triangulation of the rectangular
% region described by the vectors xvec and yvec. The entries with xvec and
% yvec define x and y partitions of the region. The arrays nxvec and nyvec
% describe the number of divisions to generate within each of the partitions.
%

[dimxvec, temp] = size(xvec) ;
[dimyvec, temp] = size(yvec) ;

x1vec = [ ] ;
y1vec = [ ] ;

% Generate x-coordinates.
for i = 1 : dimxvec-1
    tempv = linspace(xvec(i) , xvec(i+1) , nxvec(i)+1)' ;
    x1vec = [x1vec ; tempv(1:nxvec(i))] ;
end
%x1vec = [ ] ;
x1vec = [x1vec ; xvec(dimxvec)] ;
[dimx1vec, temp] = size(x1vec);

%% Generate y-coordinates.
for i = 1 : dimyvec-1
    tempv = linspace(yvec(i) , yvec(i+1) , nyvec(i)+1)' ;
    y1vec = [y1vec ; tempv(1:nyvec(i))] ;
end
%y1vec = [ ] ;
y1vec = [y1vec ; yvec(dimyvec)] ; 
[dimy1vec, temp] = size(y1vec);

%% Set up up the coordinates of the nodes.

nodeco = [ ] ;
for i = 1 : dimy1vec-1
    nodeco = [nodeco ; x1vec, y1vec(i)*ones(dimx1vec,1); ]  ;
end
nodeco = [nodeco ; x1vec, y1vec(dimy1vec)* ones(dimx1vec,1)] ;
[nnds, temp] = size(nodeco) ;

%% Now for the definition of the triangles and their neighbours.
ntri = 0 ;
for j = 1 : dimy1vec-1
    for i = 1 : dimx1vec-1
       ll = (j-1)*dimx1vec + i ;
       lr = ll + 1 ;
       ul = j*dimx1vec + i ;
       ur = ul + 1 ;

% Row type 1
       if ( rem(j,2) == 1 )
% Column type 1
       if ( rem(i,2) == 1 )
           ntri = ntri + 1;
           elnode(ntri, 1) = ll ;
           elnode(ntri, 2) = ul ;
           elnode(ntri, 3) = ur ;
%
           if ( j < dimy1vec-1 )
               elnbgh(ntri,1) = ntri + 2*(dimx1vec-1) ;
           else
               elnbgh(ntri,1) = 0 ;
           end
           elnbgh(ntri,2) = ntri + 1 ;
           if (i == 1) 
               elnbgh(ntri,3) = 0 ;
           else
               elnbgh(ntri,3) = ntri - 1 ;
           end
%
           ntri = ntri + 1;
           elnode(ntri, 1) = ll ;
           elnode(ntri, 2) = lr ;
           elnode(ntri, 3) = ur ;
%
           if ( i < dimx1vec-1 )
               elnbgh(ntri,1) = ntri + 1 ;
           else
               elnbgh(ntri,1) = 0 ;
           end
           elnbgh(ntri,2) = ntri - 1 ;
           if (j == 1) 
               elnbgh(ntri,3) = 0 ;
           else
               elnbgh(ntri,3) = ntri - 2*(dimx1vec-1) ;
           end
% Column type 2
       else
           ntri = ntri + 1;
           elnode(ntri, 1) = lr ;
           elnode(ntri, 2) = ll ;
           elnode(ntri, 3) = ul ;
%
           elnbgh(ntri,1) = ntri - 1 ;
           elnbgh(ntri,2) = ntri + 1 ;
           if (j == 1) 
               elnbgh(ntri,3) = 0 ;
           else
               elnbgh(ntri,3) = ntri - 2*(dimx1vec-1) ;
           end
%
           ntri = ntri + 1;
           elnode(ntri, 1) = lr ;
           elnode(ntri, 2) = ur ;
           elnode(ntri, 3) = ul ;
%
           if ( j < dimy1vec-1 )
               elnbgh(ntri,1) = ntri + 2*(dimx1vec-1) ;
           else
               elnbgh(ntri,1) = 0 ;
           end
           elnbgh(ntri,2) = ntri - 1 ;
           if (i < dimx1vec-1 ) 
               elnbgh(ntri,3) = ntri + 1 ;
           else
               elnbgh(ntri,3) = 0 ;
           end
       end
% Row type 2
       else
% Column type 1
       if ( rem(i,2) == 1 )
           ntri = ntri + 1;
           elnode(ntri, 1) = ul ;
           elnode(ntri, 2) = ll ;
           elnode(ntri, 3) = lr ;
%
           elnbgh(ntri,1) = ntri - 2*(dimx1vec-1) ;
           elnbgh(ntri,2) = ntri + 1 ;
           if (i == 1) 
               elnbgh(ntri,3) = 0 ;
           else
               elnbgh(ntri,3) = ntri - 1 ;
           end
%
           ntri = ntri + 1;
           elnode(ntri, 1) = ul ;
           elnode(ntri, 2) = ur ;
           elnode(ntri, 3) = lr ;
%
           if ( i < dimx1vec-1 )
               elnbgh(ntri,1) = ntri + 1 ;
           else
               elnbgh(ntri,1) = 0 ;
           end
           elnbgh(ntri,2) = ntri - 1 ;
           if (j < dimy1vec-1 ) 
               elnbgh(ntri,3) = ntri + 2*(dimx1vec-1) ;
           else
               elnbgh(ntri,3) = 0 ;
           end
% Column type 2
       else
           ntri = ntri + 1;
           elnode(ntri, 1) = ur ;
           elnode(ntri, 2) = ul ;
           elnode(ntri, 3) = ll ;
%
           elnbgh(ntri,1) = ntri - 1 ;
           elnbgh(ntri,2) = ntri + 1 ;
           if ( j < dimy1vec-1 )
               elnbgh(ntri,3) = ntri + 2*(dimx1vec-1) ;
           else
               elnbgh(ntri,3) = 0 ;
           end
%
           ntri = ntri + 1;
           elnode(ntri, 1) = ur ;
           elnode(ntri, 2) = lr ;
           elnode(ntri, 3) = ll ;
%
           elnbgh(ntri,1) = ntri - 2*(dimx1vec-1) ;
           elnbgh(ntri,2) = ntri - 1 ;
           if (i < dimx1vec-1 ) 
               elnbgh(ntri,3) = ntri + 1 ;
           else
               elnbgh(ntri,3) = 0 ;
           end
       end
       end
    end
end

%  Initially we assume that out approximation is piecewise linear
%  on all the triangles.
eldegr = ones(ntri,1) ;

%  All triangles start at level 0
ellevels = zeros(ntri,1) ;

%  Create the list that tracks the triangles currently in use
curels = zeros(ntri,1) ;
for i=1:ntri,
    curels(i) = i ;
end

% Now we define the boundary nodes and the associated boundary conditions.
% If the node has Dirichlet conditions on the intervals either side: 0
% Dirichlet condition to the left, Neumann to the right: 1
% Neumann to the left and right: 2
% Neumann to the left and Dirichlet to the right: 3
%
% In the case presented here we have that the region represents a 
% rectangular region given by:
%       [xvec(1), xvec(dimxvec)]x[yvec(1), yvec(dimyvec)],
% with Dirichlet boundary conditions all around.
%
%*** NOTE: The boundary is defined counterclockwise.
nbdy = 0 ;
for i = 1:dimx1vec
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = i ;
    bdynde(nbdy,2) = 0 ;
end
for j = 2:dimy1vec
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = j*dimx1vec ;
    bdynde(nbdy,2) = 0 ;
end
for i = 1:dimx1vec-1
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = nnds - i ;
    bdynde(nbdy,2) = 0 ;
end
for j = dimy1vec-1:-1:2
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = (j-1)*dimx1vec + 1 ;
    bdynde(nbdy,2) = 0 ;
end
