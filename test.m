clear variables;
close all;
clc;

nIt = 100;
dt = 0.1;

A = [1 dt; 0 1];
B = [dt^2/2; dt];
C = [1 0];
D = 0;
Q = [1 0; 0 1];
R = 5;
x0 = [0; 10];
x = linspace(1, 10*pi, nIt);
u = ones(nIt,1);

uscitaSigma = zeros(size(C,1), nIt);
statoSistema = zeros(length(x0), nIt);
statoKalman = zeros(length(x0), nIt);

sigma = sistema(A, B, C, D, Q, R, x0);
kal = kalmanfilter(A, B, C, D, Q , R, [0; 0]);


for i=1:nIt
    
    sigma.update(u(i));
    kal.update(u(i), sigma.leggiUscita());
        
    uscitaSigma(:,i)=sigma.leggiUscita();
    statoSistema(:,i)=sigma.leggiStato();
    statoKalman(:,i)=kal.leggiStato();
      
end

hold on;
plot(statoSistema(1, :), 'r.');
plot(uscitaSigma(1, :), 'b-');
plot(statoKalman(1, :), 'g');
legend('stato sistema','uscita sensore','stima di Kalman');
