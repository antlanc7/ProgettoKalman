classdef kalmanfilter < sistema
    %KALMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigmaModel, P, K;
    end
    
    methods
        function obj = kalmanfilter(sigmaModel, x0)
            % inizializzo la stima come l'uscita del sistema da osservare
            obj@sistema(sigmaModel.A, sigmaModel.B, sigmaModel.C, sigmaModel.Q, sigmaModel.R, x0);
            obj.sigmaModel = sigmaModel;
            obj.P = eye(obj.n);
        end
        
        function kalmanGain(obj)
            % calcolo della matrice di guadagno di kalman
            obj.K = obj.P*obj.C'*inv(obj.C*obj.P*obj.C'+obj.R);
        end
        
        function update(obj, u, y)
            %predizione
            obj.x = obj.A*obj.x + obj.B*u;
            obj.P = obj.A*obj.P*obj.A'+obj.Q;
            
            obj.kalmanGain();
            
            %aggiornamento
            obj.x = obj.x+obj.K*(y-obj.C*obj.x);
            obj.P = (eye(obj.n)-obj.K*obj.C)*obj.P;
            
        end
            
        function y = leggiUscita(obj)
            y = obj.x;
        end
        
    end
end

