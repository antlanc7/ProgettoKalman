classdef sistema < handle
    %SISTEMA 
    %Classe che descrive un sistema dinamico tempo discreto invariante
    
    properties (Access = protected)
        A,B,C,D,W,Q,R,x;    %A,B,C,D matrici del sistema
        % W matrice di guadagno del rumore di processo
        % Q matrice di covarianza del rumore di processo
        % R matrice di covarianza del rumore di misura
        n,m,p,q;      %n dim stato, m dim ingresso, p dim uscita, q dim rumore di processo
        u;       %ultimo ingresso ricevuto
    end
    
    methods
        function obj = sistema(A,B,C,D,W,Q,R,x0)
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
                obj.u = zeros(obj.m,1);
            else
                error("Matrice B non coerente con A");
            end
            
            if (size(C,2)==obj.n) % controlla che il numero di colonne di C sia n
                obj.C = C;
                obj.p = size(C,1); % p = dimensione del vettore di uscita
            else
                error("Matrice C non coerente con A");
            end
            
            if (size(D, 1) == obj.p && size(D,2) == obj.m)
                obj.D = D;
            else
                error("Matrice D non coerente con A");
            end
            
            if (size(W,1)==obj.n) % controlla che il num di righe di W sia n
                obj.W = W;
                obj.q = size(B,2); % q = dimensione del vettore rumore di processo
            else
                error("Matrice W non coerente con A");
            end
        
            
            if (isequal(size(Q),[obj.q obj.q])) % controlla che Q sia quadrata e coerente con W
                if all(abs(eig(Q))>0)
                    obj.Q = Q;
                else
                    error("Matrice Q non definita positiva");
                end
            else
                error("Matrice Q non coerente con W");
            end
            
            if (isequal(size(R),[obj.p obj.p])) % controlla che R sia quadrata e della stessa dimensione dell' uscita 
               if all(eig(R)>0)
                    obj.R = R;
               else
                   error("Matrice R non definita positiva");
               end
            else
                error("Matrice R non p*p");
            end
           
            if (isequal(size(x0),[1 obj.n])) %se x0 e' un vettore riga, delle dimensioni giuste, lo traspongo
                obj.x=x0';
            elseif (isequal(size(x0),[obj.n 1])) %se x0 e' un vettore colonna delle dimensioni giuste, OK
                obj.x=x0;
            else
                error("x0 non e' un vettore coerente con A");
            end
            
        end
                
                
        function update(obj, u)
            % aggiorna lo stato del sistema salvando quello vecchio nel vettore xold per plottare
            if (nargin<2)
                u = zeros(obj.m,1);  % se u viene omesso si considera nullo
            end
            obj.u=u;
            xn=obj.A*obj.x + obj.B*obj.u + obj.W*mvnrnd(zeros(obj.q,1),obj.Q)'; % calcola il nuovo stato x(t) = Ax(t-1) + Bu + Ww : w = rumore di processo
            obj.x = xn;               % aggiorna lo stato con quello nuovo
            
        end
        
        function y = leggiUscita(obj)
            % restituisce in output l'uscita y del sistema
            y=obj.C*obj.x + obj.D*obj.u + mvnrnd(zeros(obj.p,1),obj.R)'; % y = Cx + v : v = rumore di misura
        end
        
        %get dello stato per plot
        function esc = leggiStato(obj)
            esc = obj.x;
        end
        
    end
end

