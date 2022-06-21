%% TDMA
function [phi,iter,ResImax,tempo_tdma] = tdma2d(ap,an,as,aw,ae,Su)
tic
[row, col] = size(ap);
% Coeficientes do TDMA
Aj = zeros(row,col);
Cj = zeros(row,col);
Clj= zeros(row,col);

% Incialização do campo da propriedade phi
phi = zeros(row, col);

% Iniciando os residuos e iterações
ResG     = [];                  % Global
ResI     = zeros(row,col);      % Matriz de residuo local
iter     = 0; 
ResImax  = []; 
% Tolerância do TDMA
tol = 1e-3;
CP  = 1; % Critério de Parada 

while(CP >= tol)
        %%%% Varredura linha-a-linha da malha sentido: Norte->Sul %%%%
    for i=1:row        
        if i == 1                   % Extremo Norte da Malha
            phi_n   = 0;
        else 
            phi_n = phi(i - 1,:);    
        end
        
        if i == row                 % Extremo Sul da Malha
            phi_s = 0;
        else
            phi_s = phi(i + 1,:);
        end   
        
        Cj(i,:) =  Su(i,:) + an(i,:).*phi_n + as(i,:).*phi_s;
        
        %%%% Varredura dos vc's na linha i no sentido: Oeste->Leste %%%%
        for j = 1:col
            k = col - (j-1);        % índice para o back substitution 
            
            if j == 1               % Extremo Oeste da Malha
                phi_w   = 0;
                Ajlast  = 0;
                Cljlast = 0;
            else                   
                phi_w   = phi(i,j-1);
                Cljlast = Clj(i,j-1);
                Ajlast  = Aj(i,j-1);
            end
            
            if j == col             % Extremo Leste da Malha
                phi_e = 0;
            else
                phi_e = phi(i,j+1);
            end
            
                    %%%% foward elimination %%%%
            Aj(i,j)  = ae(i,j)/(ap(i,j) - aw(i,j)* Ajlast);
            Clj(i,j) = (Cj(i,j) + aw(i,j)*Cljlast)/(ap(i,j) - aw(i,j)*Ajlast);
        
                    %%%% backward substitution %%%%
            if k == col
                phi(i,k) = Clj(i,k);
            else
                phi(i,k) = Aj(i,k)*phi(i,k+1) + Clj(i,k);
            end
            
                    %%%% Calculo os Residuos Locais %%%%
            aux         = aw(i,j)*phi_w + ae(i,j)*phi_e;
            ResI(i,j)   = Cj(i,j) + aux - ap(i,j)*phi(i,j);    
        end 
    end
    
                %%%% Atualiza o Critério de Parada %%%%
    CP  = max(max(ResI.^2)); 
    
        %%%% Atualiza a iteração e Calculo dos Residuos %%%%
    iter            = iter + 1; 
    ResG(iter)      = sum(sum(abs(ResI)));
    ResImax(iter)   = CP;
end
tempo_tdma = toc;

end