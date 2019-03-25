classdef sistema
    %SISTEMA 
    %Classe che descrive un sistema dinamico tempo discreto invariante
    
    properties
        A,B,C,Q,R,x;    %A,B,C matrici del sistema
        % Q matrice di covarianza del rumore di processo
        % R matrice di covarianza del rumore di misura
        n,m,p;      %n dim stato, m dim ingresso, p dim uscita
    end
    
    methods
        function obj = sistema(A,B,C,Q,R,x0)
            %SISTEMA Construct an instance of this class
            %   Detailed explanation goes here
            if (diff(size(A))==0) 
                obj.A = A;
                obj.n=length(A);
            else
                error("ERRORE: Matrice A non quadrata");
                return
            end
            obj.B = B;
            obj.C = C;
            obj.Q = Q;
            obj.R = R;
            
            obj.x(1)=x0;
            
        end
        
        function y = leggiUscita(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            y=obj.C*obj.x;
        end
    end
end

