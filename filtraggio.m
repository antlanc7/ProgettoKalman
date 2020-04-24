close all; clear; clc;

dt=1e-2;        % passo di campionamento considerato
durata=10;      % durata simulazione
t=0:dt:durata;  % costruzione vettore dei tempi

% dimensioni sistema
n=3;
m=1;
p=1;
A=zeros(n);
for i=1:n-1
    A(i,i+1)=1;
end
A
%A=[0 1; -4*pi*pi 0];
B=zeros(n,m);
B(n)=1
C=zeros(p,n);
C(1)=1
D=zeros(p,m);
R=1e-2*eye(p);

%discretizzazione
sys=ss(A,B,C,D);
sysd=c2d(sys,dt);
[Ad,Bd,Cd,Dd]=ssdata(sysd)

Q=1e-5

opt=menu("Scegli segnale da filtrare","scalino","rampa","parabola","esponenziale convergente","sinusoide","sinusoide smorzata","onda quadra");
switch (opt)
    case 0, return
    case 1, x = ones(1,length(t));
    case 2, x = t/durata;
    case 3, x = t.^2/durata^2;
    case 4, x = exp(-t/2);
    case 5, x = 1*sin(2*pi*t);
    case 6, x = exp(-t/2).*sin(2*pi*t);
    case 7, x = square(2*pi*t);
end

%generazione rumore bianco
v=mvnrnd(0,R,length(t))';

% aggiunta rumore al segnale
xe=x+v;

% inizializzazione filtro di kalman
x0=zeros(n,1);
P0=eye(n);
k=kalmanfilter(Ad,Bd,Cd,Dd,Q,R,x0,P0);

% salvataggio valori matrici K e P per plot
Kplot=zeros(n,p,length(t));
Pplot=zeros(n,n,length(t));

xf=zeros(1,length(t));
%xpr=zeros(1,length(t)); vettore delle x predette

for i=1:length(t)
    xf(i)=k.leggiUscita();
    %xpr(i)=k.leggiXpr();
    k.update(0,xe(i));
    Kplot(:,:,i)=k.leggiK();
    Pplot(:,:,i)=k.leggiP();
end

figure('WindowState','maximized');
subplot(2,3,[1 2 3]);
title("Filtro di Kalman");
hold on;
grid on;
plot(t,xe,'b.') % campioni
plot(t,xf(1,:),'r','LineWidth',2); % stima kalman
plot(t,x,'k') % segnale puro
%plot(t,xf(1,:)-x); % errore
legend("Segnale campionato con rumore","Stima del filtro di Kalman","Segnale originale");

subplot(2,3,4);
title("Tempo t in secondi");
grid on;
tline=animatedline;
xlim([0 durata]);
ylim([-1 1]);

for i=1:durata/dt/100:length(t)
    subplot(2,3,4);
    addpoints(tline,t(i),0);
    
    subplot(2,3,5);
    bar3(Kplot(:,:,i));
    title("Matrice K");
    
    subplot(2,3,6);
    bar3(Pplot(:,:,i));
    title("Matrice P");
    
    pause(dt)
end
