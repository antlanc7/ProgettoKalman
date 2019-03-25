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
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.Q = Q;
            obj.R = R;
            
            obj.x(1)=x0;
            
        end
        
        function y = leggiUscita()
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            y=C*x+randn();
        end
        
        function update(obj, u)
          obj.x(end+1) = obj.A*obj.x(end) + obj.B*u;
        end
        
    end
end

