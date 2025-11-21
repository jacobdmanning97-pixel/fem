
%% To Plot the x and y components of the velocity.
TEl = [elnode(:,1:3) , nVert + elnode(:,4:6)] ;
nplt = size(GlobalV,1) ; 

if ( strcmp(vel_bas_type, 'CtsQuad') == 1 ) 
   	trisurf([TEl(:,[1 6 5]); TEl(:,[2 4 6]); TEl(:,[3 5 4]); TEl(:,[4 5 6])], ...
   		nodeco(:,1), nodeco(:,2), GlobalV(1:nplt,1)), shading interp,  ...
   		colorbar

elseif strcmp(vel_bas_type, 'CtsLin') == 1
   	trisurf(TEl(:,1:3), ...
   		nodeco(1:nVert,1), nodeco(1:nVert,2), GlobalV(1:nplt,1)), shading interp,  ...
   		colorbar
end
   


