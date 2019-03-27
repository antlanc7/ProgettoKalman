clear all
close all
clc

dt = 0.1;

A = [1 dt; 0 1];
B = [dt*dt/2; dt];
C = [1 0];
Q = [0 0; 0 1];
R = 7;
x0 = [50; 10];
u = 1;

sigma = sistema(A, B, C, Q, R, x0);
kal = kalmanfilter(sigma, [0; 0]);

for i=1:100
        
    uscitaSigma(:,i) = sigma.leggiUscita();
    statoSistema(:,i)=sigma.leggiStato();
    statoKalman(:,i)=kal.leggiStato();
  
    sigma.update(u);
    kal.update(u, sigma.leggiUscita());
    
end

hold on;
plot(statoSistema(1, :), 'r.');
plot(statoKalman(1, :), 'g');
plot(uscitaSigma(1, :), 'b-');
legend('stato sistema','Stima di Kalman','uscita sistema');
