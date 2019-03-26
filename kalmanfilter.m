classdef kalmanfilter < handle
    %KALMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigmaModel, P, K, stima;
    end
    
    methods
        function obj = kalman(sigmaModel)
            % inizializzo la stima come l'uscita del sistema da osservare
            obj.sigmaModel = sigmaModel;
            obj.P = inv(sigmaModel.C)*sigmaModel.R*inv(sigmaModel.C');
            obj.stima = sigmaModel.leggiUscita();
        end
        
        function kalmanGain(obj)
            % calcolo della matrice di guadagno di kalman
            obj.K = obj.P*obj.sigmaModel.C'*inv(obj.sigmaModel.C*obj.P*obj.sigmaModel.C'+obj.sigmaModel.R);
        end
        
    end
end

