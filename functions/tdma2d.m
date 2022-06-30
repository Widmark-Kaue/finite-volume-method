%% TDMA (Tri Diagonal Matrix Algorithm)
%{
INSTRUCOES:
------------
Para utilizar essa função no MATLAB do terminal remoto da USFC este
arquivo de função <tdma2d.m> deve ser copiado para o seguinte 
diretório:

/Este Computador/Documentos/MATLAB/

%}
function [phi,iter, Res, tempo_tdma] = tdma2d(ap,an,as,aw,ae,Su)
tic
[row, col] = size(ap);

% Coeficientes do TDMA
Aj = zeros(row,col);
Cj = zeros(row,col);
Clj= zeros(row,col);

% Incializacao do campo da propriedade phi
phi = zeros(row, col);

% Iniciando os residuos e iteracoes
ResG     = [];                  % Global
ResI     = zeros(row,col);      % Matriz de residuo local
iter     = 0; 
ResImax  = []; 
ResN     = [];
% Tolerancia do TDMA
mi_ap = min(min(ap));
if mi_ap < 1
    tol = 1e-3*mi_ap;
else
    tol = 1e-3;
end
CP1  = 1; % Criterio de Parada 
CP2  = 1;
while(CP1 >= tol || CP2 >= tol)
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
    
    
        %%%% Atualiza a iteracao e Calculo dos Residuos %%%%
    iter            = iter + 1; 
    ResG(iter)      = sum(sum(abs(ResI)));
    
                %%%% Atualiza o Criterio de Parada %%%%
    CP1  = max(max(ResI.^2)); 
    if iter >= 5
        CP2 = ResG(iter)/ResG(5);
    end
    ResImax(iter)   = CP1;
    ResN(iter)      = CP2;
end
Res         =[ResImax; ResN];
tempo_tdma  = toc;
%adicionar método de parada por número de iteraçao
end
