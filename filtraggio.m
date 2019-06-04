close all; clearvars; clc;

dt=1e-3;
durata=1.5;
t=0:dt:durata;

n=3;
A=zeros(n);
for i=1:n-1
    A(i,i+1)=1;
end
A
B=zeros(n,1);
B(n)=1
C=zeros(1,n);
C(1)=1
D=0;
R=1e-2;

eAdt=expm(A*dt)
Q=R*integral(@(x)expm(A*x)*(B*B')*expm(A'*x),0,dt,'ArrayValued',true)

opt=menu("Scegli segnale da filtrare","scalino","rampa","parabola","sinusoide","esponenziale convergente","sinusoide smorzata","onda quadra");
switch (opt)
    case 0, return
    case 1, x = t.^0;
    case 2, x = t;
    case 3, x = t.^2/2;
    case 4, x = sin(200*pi*t);
    case 5, x = exp(-t/2);
    case 6, x = exp(-t/2).*sin(2*pi*t);
    case 7, x = (1 + sign(cos(pi*t)))/2;
end

xe=zeros(1,length(t));
% aggiunta rumore
for i=1:length(t)
    xe(i) = x(i) + mvnrnd(0,R);
end

k=kalmanfilter(eAdt,B,C,D,Q,R,zeros(n,1),eye(n));

xf=zeros(n,length(t));
for i=1:length(t)
    xf(:,i)=k.leggiUscita();
    k.update(0,xe(i));
end

hold on
plot(t,x,'k')
plot(t,xe,'b.')
plot(t,xf(1,:),'r','LineWidth',2);
plot(t,xf(1,:)-x);
