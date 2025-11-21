%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type
global quad_rul num mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vel_bas_type = 'CtsLin' ;
quad_rul = 'quad_73';
num = 1;
listmesh = ["R4x4","R6x6","R8x8","R10x10","R12x12"];

l = zeros(5,5);
for i=1:5
    mesh = listmesh(i);
    [uL2Error, GraduL2Error, H1uError]=Drv_ConDiff;
    l(i,:)=[uL2Error,2*i+2, 0,H1uError, 0];
end
for i=1:4
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i+1,2)/l(i,2));
    l(i+1,5) = log(l(i,4)/l(i+1,4))/log(l(i+1,2)/l(i,2));
end
lin=table(l(:,2),l(:,1),l(:,3),l(:,4),l(:,5));
lin.Properties.VariableNames=["Mesh","uL2Error","alpha uL2","H1uError","alpha H1u"];

disp("Part 2")
disp("Using Continuous Linears")
disp(lin)
disp("Thus the experimental convergence rate is about 1 for H_1 norm, which is what we would expect since" +...
    " the error in the H_1 norm should be in the first term because we are using continuous linear aproximations" + ...
    " and the error for the L2 norm is about 2 which is again what we would expect using continuous linears.")
disp("I believe that my program works correctly as it solves the given problems to the error" + ...
    " that we would expect given the continuous linear approximation.")

vel_bas_type = 'CtsQuad' ;
l = zeros(5,5);
for i=1:5
    mesh = listmesh(i);
    [uL2Error, GraduL2Error, H1uError]=Drv_ConDiff;
    l(i,:)=[uL2Error,2*i+2, 0,H1uError, 0];
end
for i=1:4
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i+1,2)/l(i,2));
    l(i+1,5) = log(l(i,4)/l(i+1,4))/log(l(i+1,2)/l(i,2));
end
quad=table(l(:,2),l(:,1),l(:,3),l(:,4),l(:,5));
quad.Properties.VariableNames=["Mesh","uL2Error","alpha uL2","H1uError","alpha H1u"];


disp("Using Continuous Quadratics")
disp(quad)

disp("I believe that my program works correctly as it solves the given problems to machine error" + ...
    " which is what we would expect given the continuous quadratic approximation on a quadratic solution.")


vel_bas_type = 'CtsLin' ;
num = 0;

l = zeros(5,5);
for i=1:5
    mesh = listmesh(i);
    [uL2Error, GraduL2Error, H1uError]=Drv_ConDiff;
    l(i,:)=[uL2Error,2*i+2, 0,H1uError, 0];
end
for i=1:4
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i+1,2)/l(i,2));
    l(i+1,5) = log(l(i,4)/l(i+1,4))/log(l(i+1,2)/l(i,2));
end
lin=table(l(:,2),l(:,1),l(:,3),l(:,4),l(:,5));
lin.Properties.VariableNames=["Mesh","uL2Error","alpha uL2","H1uError","alpha H1u"];

disp("Part 3a")
disp("Using Continuous Linears")
disp(lin)
disp("Thus the experimental convergence rate is about 1 for H_1 norm, which is what we would expect since" +...
    " the error in the H_1 norm should be in the first term because we are using continuous linear aproximations" + ...
    " and the error for the L2 norm is about 2 which is again what we would expect using continuous linears.")

vel_bas_type = 'CtsQuad' ;
l = zeros(5,5);
for i=1:5
    mesh = listmesh(i);
    [uL2Error, GraduL2Error, H1uError]=Drv_ConDiff;
    l(i,:)=[uL2Error,2*i+2, 0,H1uError, 0];
end
for i=1:4
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i+1,2)/l(i,2));
    l(i+1,5) = log(l(i,4)/l(i+1,4))/log(l(i+1,2)/l(i,2));
end
quad=table(l(:,2),l(:,1),l(:,3),l(:,4),l(:,5));
quad.Properties.VariableNames=["Mesh","uL2Error","alpha uL2","H1uError","alpha H1u"];


disp("Using Continuous Quadratics")
disp(quad)

disp("Thus the experimental convergence rate is about 2 for H_1 norm, which is what we would expect since" +...
    " the error in the H_1 norm should be in the second term because we are using continuous quadratic aproximations" + ...
    " and the error for the L2 norm is about 3 which is again what we would expect using continuous quadratics.")


quad_rul = 'quad_75';

l = zeros(5,5);
for i=1:5
    mesh = listmesh(i);
    [uL2Error, GraduL2Error, H1uError]=Drv_ConDiff;
    l(i,:)=[uL2Error,2*i+2, 0,H1uError, 0];
end
for i=1:4
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i+1,2)/l(i,2));
    l(i+1,5) = log(l(i,4)/l(i+1,4))/log(l(i+1,2)/l(i,2));
end
quad=table(l(:,2),l(:,1),l(:,3),l(:,4),l(:,5));
quad.Properties.VariableNames=["Mesh","uL2Error","alpha uL2","H1uError","alpha H1u"];

disp("Part 3b")
disp("Using Continuous Quadratics with quad_75")
disp(quad)

disp("I did not see a difference in the experimental convergence rate after changing the quad rule." + ...
    " If anything, it got worse.")



