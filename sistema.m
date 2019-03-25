classdef sistema < handle
    %SISTEMA 
    %Classe che descrive un sistema dinamico tempo discreto invariante
    
    properties
        A,B,C,Q,R,x;    %A,B,C matrici del sistema
        % Q matrice di covarianza del rumore di processo
        % R matrice di covarianza del rumore di misura
        n,m,p;      %n dim stato, m dim ingresso, p dim uscita
        xold;       %vettore stati vecchi (per plot)
    end
    
    methods
        function obj = sistema(A,B,C,Q,R,x0)
            %SISTEMA Costruttore
            
            if (diff(size(A))==0) % se A è quadrata
                obj.A = A;
                obj.n=size(A,1);
            else
                error("Matrice A non quadrata");
                return
            end
            
            if (size(B,1)==obj.n)
                obj.B = B;
                obj.m = size(B,2);
            else
                error("Matrice B non coerente con A");
                return
            end
            
            if (size(C,2)==obj.n)
                obj.C = C;
                obj.p = size(C,1);
            else
                error("Matrice C non coerente con A");
                return
            end
            
            if (diff(size(Q))==0 && size(Q,1)==obj.n) 
                obj.Q = Q;
            else
                error("Matrice Q non coerente con x");
                return
            end
            
            if (diff(size(R))==0 && size(R,1)==obj.p) 
                obj.R = R;
            else
                error("Matrice R non coerente con y");
                return
            end
           
            if (isequal(size(x0),[1 obj.n]))
                obj.x=x0';
            elseif (isequal(size(x0),[obj.n 1]))
                obj.x=x0;
            else
                error("Vettore x0 non coerente con A");
                return
            end
        end
        
        function y = leggiUscita(obj)
            % restituisce in output l'uscita y del sistema
            y=obj.C*obj.x + obj.R*randn(obj.p,1);
        end
        
        function update(obj, u)
            % aggiorna lo stato del sistema salvando quello vecchio nel
            % vettore xold per plottare
            obj.xold(:,end+1)=obj.x;  % salva il vecchio stato
            xn=obj.A*obj.x + obj.B*u + obj.Q*randn(obj.n,1); % calcola il nuovo stato 
            obj.x = xn;               % aggiorna lo stato con quello nuovo
        end
        
    end
end

