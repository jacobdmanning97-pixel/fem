format long
type = 'd1_CtsLin' ;

alpha = 0.0;
beta = 2.0;

[uErr, DuErr] = Driver_1d(11,alpha,beta,true,1,type);
disp("Part A 2bi")
disp("Given the solution in the test space u(x)=3x-3 and 11 nodes, uErr= "+ uErr +" and DuErr= " +DuErr)

nnds = 2.^(3:7)+1 ;
l = zeros(size(nnds,2),3);
for i=1:size(nnds,2)
    h = (beta - alpha) / (nnds(i) - 1) ;
    %Driver_1d(nnds,alpha,beta, uTrue_Flag,num)
    [uErr, DuErr] = Driver_1d(nnds(i),alpha,beta,true,3,type);
    l(i,:)=[DuErr, h, 0];
end
for i=1:size(nnds, 2)-1
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i,2)/l(i+1,2));
end
t=table(l(:,1),l(:,2),l(:,3));
t.Properties.VariableNames=["DuErr","h","alpha DuErr"];
disp("Part A 2bii")
disp("The experimental convergence rate for x^4+1 is given in the last column of the table")
disp(t)
disp("Thus the experimental convergence rate is about 1, which is what we would expect since" +...
    " the error in the H_1 norm should be in the first term because we are using continuous linear aproximations.")
disp("I believe that my program works correctly as it solves the given problems to the error" + ...
    " that we would expect given the continuous linear approximation")

type = 'd1_CtsQuad' ;
[uErr, DuErr] = Driver_1d(11,alpha,beta,true,2,type);
disp("Part B 2bi")
disp("Given the solution in the test space u(x)=3x^2-3 and 11 nodes, uErr= "+ uErr +" and DuErr= " +DuErr)

l = zeros(size(nnds,2),3);
for i=1:size(nnds,2)
    h = (beta - alpha) / (nnds(i) - 1) ;
    %Driver_1d(nnds,alpha,beta, uTrue_Flag,num)
    [uErr, DuErr] = Driver_1d(nnds(i),alpha,beta,true,3,type);
    l(i,:)=[DuErr, h, 0];
end
for i=1:size(nnds, 2)-1
    l(i+1,3) = log(l(i,1)/l(i+1,1))/log(l(i,2)/l(i+1,2));
end
t=table(l(:,1),l(:,2),l(:,3));
t.Properties.VariableNames=["DuErr","h","alpha DuErr"];
disp("Part B 2bii")
disp("The experimental convergence rate for x^4+1 is given in the last column of the table")
disp(t)
disp("Thus the experimental convergence rate is about 2, which is what we would expect since" +...
    " the error in the H_1 norm should be in the second term because we are using continuous quadratic aproximations.")
disp("I believe that my program works correctly as it solves the given problems to the error" + ...
    " that we would expect given the continuous quadratic approximation")
