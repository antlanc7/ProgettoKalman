classdef kalmanfilter < sistema
    
    %classe filtro di Kalman:
    %eredita da sistema essendo a sua volta un sistema
    
    properties (Access = protected)
        P, K;   %P matrice di covarianza dello stato; K matrice guadagno di Kalman
    end
    
    methods
        
        %costruttore della classe: inizializza il filtro
        function obj = kalmanfilter(sigmaModel, x0) %parametri: sistema da osservare, "stima" iniziale dello stato
            
            %chiama il costruttore della superclasse (sistema) copiando le
            %matrici del sistema da osservare (sigmaModel)
            obj@sistema(sigmaModel.A, sigmaModel.B, sigmaModel.C, sigmaModel.D, sigmaModel.Q, sigmaModel.R, x0);
            
            %inizializza la matrice di covarianza dello stato 
            %(matrice quadrata della dimensione dello stato)
            obj.P = eye(obj.n);
            
        end
        
        %calcola la stima dello stato del sistema osservato
        function update(obj, u, y) %parametri: ingresso dato a sigma(u), uscita del sistema osservato(y)
            
            %predizione
            obj.x = obj.A*obj.x + obj.B*u;
            obj.P = obj.A*obj.P*obj.A'+obj.Q;
            
            %calcolo del guadagno di Kalman
            obj.K = obj.P*obj.C'/(obj.C*obj.P*obj.C'+obj.R);
          
            %correzione
            obj.x = obj.x+obj.K*(y-obj.C*obj.x);
            obj.P = (eye(obj.n)-obj.K*obj.C)*obj.P;
            
        end
        
        function y = leggiUscita(obj)
            y = obj.x;  %uscita dell'osservatore uguale allo stato
        end
        
    end
end

