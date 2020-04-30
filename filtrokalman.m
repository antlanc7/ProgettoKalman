classdef filtrokalman < sistema
    
    %classe filtro di Kalman:
    %eredita da sistema essendo a sua volta un sistema
    
    properties %(Access = protected)
        L;          % matrice guadagno di Kalman
        P;          % matrice di covarianza dello stato corretto
        xPr, PPr;   % predizione dello stato e relativa covarianza
    end
    
    methods
        
        %costruttore della classe: inizializza il filtro
        function obj = filtrokalman(A, B, C, D, W, Q, R, x0, P0) %parametri: sistema da osservare, "stima" iniziale dello stato
            
            %chiama il costruttore della superclasse (sistema) passando le
            %matrici
            obj@sistema(A, B, C, D, W, Q, R, x0);
            
            %inizializza la matrice di covarianza dello stato 
            %(matrice quadrata della dimensione dello stato)
            if (nargin<8)
                P0 = eye(obj.n);  % se P0 viene omessa si considera l'identità
            end
            
            obj.xPr=x0;
            obj.PPr=P0;
        end
        
        %calcola la stima dello stato del sistema osservato
        function update(obj, u, y) %parametri: ingresso (u) e uscita (y) del sistema osservato
            
            obj.u=u;
            
            %calcolo guadagno di Kalman
            obj.L = obj.PPr*obj.C'/(obj.C*obj.PPr*obj.C'+obj.R);
          
            %correzione
            obj.x = obj.xPr+obj.L*(y-obj.C*obj.xPr);
            I_LC = (eye(obj.n)-obj.L*obj.C);
            obj.P = I_LC*obj.PPr*I_LC'+obj.L*obj.R*obj.L';
            
            %predizione
            obj.xPr = obj.A*obj.x + obj.B*u;
            obj.PPr = obj.A*obj.P*obj.A'+obj.W*obj.Q*obj.W';
            
        end
        
        function x = leggiStima(obj)
            x = obj.x;
        end
        
        function y = leggiUscitaStimata(obj)
            y = obj.C*obj.x + obj.D*obj.u;
        end
        
        function L = leggiL(obj)
            L = obj.L;
        end
        
        function P = leggiP(obj)
            P = obj.P;
        end
        
    end
end

