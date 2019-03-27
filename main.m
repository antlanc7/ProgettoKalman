clear all
close all
clc

dt = 0.1;

A = [1 dt; 0 1];
B = [dt*dt/2; dt];
C = [1 0];
Q = [0 0; 0 1];
R = 1;
x0 = [0; 0];
u = 1;

sigma = sistema(A, B, C, Q, R, x0);
kal = kalmanfilter(sigma, [0; 0]);

for i=1:100
   % u = i;
    sigma.update(u);
    y = sigma.leggiUscita();
    kal.update(u, y);
    
    uscitaSigma(:,i) = y;
    statoSistema(:,i)=sigma.x;
    statoKalman(:,i)=kal.x;
end

hold on;
plot(statoSistema(1, :), 'r.');
plot(statoKalman(1, :), 'g');
plot(uscitaSigma(1, :), 'b--');
legend('stato sistema','Kalman output','uscita sistema');
