%Problem 2b
e2=1/20;

x21=fd(e2,20);
x22=fd(e2,40);

figure
hold on
fplot(@(x) exp(1/e2)/(exp(1/e2)-1)-exp(x/e2)/(exp(1/e2)-1))
plot(linspace(0,1,21), x21)
plot(linspace(0,1,41), x22)
axis([0 1 0 2])
hold off

legend('True Solution','N=20','N=40')

%Problem 5
e5=1/200;

x51=fd(e5,20);
x52=fd(e5,40);

figure
hold on
fplot(@(x) exp(1/e5)/(exp(1/e5)-1)-exp(x/e5)/(exp(1/e5)-1))
plot(linspace(0,1,21), x51)
plot(linspace(0,1,41), x52)
axis([0 1 0 2])
hold off

legend('True Solution','N=20','N=40')

%Problem 6

x61=fd(e5,1000);
x62=fe(e5,1000);

figure
hold on
fplot(@(x) exp(1/e5)/(exp(1/e5)-1)-exp(x/e5)/(exp(1/e5)-1))
plot(linspace(0,1,1001), x61)
plot(linspace(0,1,1001), x62)
axis([0 1 0 2])
hold off

legend('True Solution','Finite Difference','Finite Element')

function xs=fd(ep,N)
    h=1/N;
    B=diag([1 2*ep/h^2*ones(1,N-1) 1])+diag((-ep/h^2+1/(2*h))*ones(1,N),1)+diag((-ep/h^2-1/(2*h))*ones(1,N),-1);
    B(1,2)=0;
    B(N+1,N)=0;
    s=zeros(N+1,1);
    s(1)=1;
    
    xs=B\s;
end

function xs=fe(ep,N)
    h=1/N;
    B=diag([1 2*ep/h*ones(1,N-1) 1])+diag((-ep/h+1/2)*ones(1,N),1)+diag((-ep/h-1/2)*ones(1,N),-1);
    B(1,2)=0;
    B(N+1,N)=0;
    s=zeros(N+1,1);
    s(1)=1;
    
    xs=B\s;
end