close all; clear; clc;

dt=1e-2        % passo di campionamento
durata=10      % durata simulazione
t=0:dt:durata; % costruzione vettore dei tempi


% menu scelta segnale
opt=menu("Selezione segnale da filtrare:","scalino","rampa","parabola","esponenziale convergente","sinusoide","sinusoide smorzata");


%% costruzione generatore di segnale e stato iniziale
if opt==0
    return
elseif opt<=3   %polinomi
    n=opt;  %dimensione sistema
    k=n-1;  %grado polinomio
    A=[zeros(k,1) eye(k); 0 zeros(1,k)]
    if opt==1
        x0=1;
    else
        x0(n)=1/(durata^k);
    end
elseif opt==4   % esponenziale
    n=1;
    alfa=-1/2;
    A=alfa
    x0=exp(1);
    %x = exp(alfa*t)
elseif opt==5   % sinusoide
    n=2;
    omega=2*pi;
    A=[0 -omega; omega 0]
    x0=[1;0];
    %x1 = cos(omega*t)
    %x2 = sin(omega*t)
elseif opt==6   % sinusoide smorzata
    n=2;
    alfa=-1/2;
    omega=2*pi;
    A = [alfa -omega; omega alfa]
    x0=[1;0];
    %x1 = cos(omega*t)*exp(alfa*t)
    %x2 = sin(omega*t)*exp(alfa*t)
else
    return
end
m=1;
p=1;

B=zeros(n,m);
C=eye(p,n);
D=zeros(p,m);

%discretizzazione matrici sistema
sys=ss(A,B,C,D);
sysd=c2d(sys,dt);
[Ad,Bd,Cd,Dd]=ssdata(sysd);

% matrici dei rumori
W=zeros(n,1);
W(n)=dt;
Q=1e-3;
R=1e-1*eye(p);

% inizializzazione sistema generatore di segnale rumoroso
sys=sistema(Ad,Bd,Cd,Dd,W,Q,R,x0);

% inizializzazione filtro di kalman
P0=eye(n);
k=filtrokalman(Ad,Bd,Cd,Dd,W,Q,R,zeros(n,1),P0);

% creazione array di matrici L e P per plot
Lplot=zeros(n,p,length(t));
Pplot=zeros(n,n,length(t));

% vettori dei segnali
x=zeros(n,length(t));
y=zeros(p,length(t));
xs=zeros(n,length(t));

%% simulazione
for i=1:length(t)
    x(:,i)=sys.leggiStato();     % lettura stato sistema (segnale da ricostruire, per plot)
    y(:,i)=sys.leggiUscita();    % lettura uscita sistema (segnale rumoroso)
    sys.update();                % calcolo del nuovo stato del sistema
    k.update(0,y(:,i));          % aggiornamento filtro --> parametri u=0 e y(segnale rumoroso)
    xs(:,i)=k.leggiStima();      % lettura stima kalman
    Lplot(:,:,i)=k.leggiL();     % lettura matrice L per animazione
    Pplot(:,:,i)=k.leggiP();     % lettura matrice P per animazione
end

%% grafica

f1=figure('WindowState','maximized');
subplot(2,2,[1 2]);
title("Filtro di Kalman");
hold on;
grid on;
plot(t,y,'b.') % campioni da filtrare
plot(t,xs(1,:),'r','LineWidth',2); % stima kalman
plot(t,x(1,:),'k') % stato sistema da ricostruire (non visibile al filtro)
%plot(t,xf(1,:)-x); % errore
legend("Segnale rumoroso","Segnale filtrato","Segnale originale");
%%
for i=1:durata/dt/100:length(t)
    if ~ishghandle(f1)
        break   %ferma l'animazione se viene chiusa la finestra
    end
    tempo=sprintf('Istante t = %.1f s', t(i));
    subplot(2,2,3);
    bar3(Lplot(:,:,i));
    title("Matrice L");
    xlabel(tempo);
    
    subplot(2,2,4);
    bar3(Pplot(:,:,i));
    title("Matrice P");
    xlabel(tempo);
    
    pause(0.03)
end


