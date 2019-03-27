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
            
            if (diff(size(A))==0) % se A e' quadrata
                obj.A = A;
                obj.n=size(A,1); % n = dimensione dello stato (e della matrice a)
            else
                error("Matrice A non quadrata");
            end
            
            if (size(B,1)==obj.n) % controlla che il num di righe di B sia n
                obj.B = B;
                obj.m = size(B,2); % m = dimensione del vettore di ingresso
            else
                error("Matrice B non coerente con A");
            end
            
            if (size(C,2)==obj.n) % controlla che il numero di colonne di C sia n
                obj.C = C;
                obj.p = size(C,1); % p = dimensione del vettore di uscita
            else
                error("Matrice C non coerente con A");
            end
            
            if (diff(size(Q))==0 && size(Q,1)==obj.n) % controlla che Q sia quadrata e della stessa dimensione dello stato
                obj.Q = Q;
            else
                error("Matrice Q non coerente con x");
            end
            
            if (diff(size(R))==0 && size(R,1)==obj.p) % controlla che R sia quadrata e della stessa dimensione dell' uscita 
                obj.R = R;
            else
                error("Matrice R non coerente con y");
            end
           
            if (isequal(size(x0),[1 obj.n])) %se x0 e' un vettore riga, delle dimensioni giuste, lo traspongo
                obj.x=x0';
            elseif (isequal(size(x0),[obj.n 1])) %se x0 e' un vettore colonna delle dimensioni giuste, OK
                obj.x=x0;
            else
                error("x0 non e' un vettore delle dimensioni giuste");
            end
        end
        
        function y = leggiUscita(obj)
            % restituisce in output l'uscita y del sistema
            y=obj.C*obj.x + obj.R*randn(obj.p,1); % y = Cx + w : w = rumore di misura
        end
        
        function update(obj, u)
            % aggiorna lo stato del sistema salvando quello vecchio nel vettore xold per plottare
            if (nargin<2) u = zeros(obj.m,1); end    % se u viene omesso si considera nullo
            obj.xold(:,end+1)=obj.x;  % salva il vecchio stato
            xn=obj.A*obj.x + obj.B*u + obj.Q*zeros(obj.n,1); % calcola il nuovo stato x(t) = Ax(t-1) + Bu + v : v = rumore di processo
            obj.x = xn;               % aggiorna lo stato con quello nuovo
        end       
    end
end

