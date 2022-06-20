clear all
close all
clc
%% Parametros do Problema
% Geometria
W   = 0.5; 
L   = 0.5;

% Propriedades do Material
k   = 1;

% Condições de contorno
T1  = 0;
T2  = 100; 
 
%% Definindo a malha

% nº de VC's nas direções x e y
nx = 50; 
ny = 50; 

% nº total de VC's
vc = nx*ny;

% Dimensoes dos VC's
dx = L/nx;
dy = W/ny;

%% Construçoẽs do Domínio (coeficientes)

row = ny; col = nx;

% Funções auxiliares
ij = @(k,n) 1 + k*(n - 1);
i  = @(k) ij(k,row); 
j  = @(k) ij(k,col); 

% Matrizes de coeficientes
Ans = k*dx/dy;
Awe = k*dy/dx;

aw = Awe*ones(row,col); 
ae = Awe*ones(row,col);
as = Ans*ones(row,col);
an = Ans*ones(row,col);
Su = zeros(row, col);


for k =0:1   
   if k == 0                % Volumes a norte e a oeste
      an(i(k),:)     = 0;
      aw(:,j(k))     = 0;
   else                     % Volumes a sul e a leste
        as(i(k),:)   = 0; 
        ae(:,j(k))   = 0;
   end
end

ap = aw + ae + as + an;

for k = 0:1
    ap(i(k),:)   = ap(i(k),:) + 2*Ans;  % Adicionando o termo Sp
    ap(:,j(k))   = ap(:,j(k)) + 2*Awe;
    
    if k == 0                           % Construção do termo fonte
        Su(i(k),:)   = Su(i(k),:) + 2*T2; 
    end
    Su(:,j(k))   = Su(:,j(k)) + 2*T1;
end
%% Solução pela inversa
tic
% Matrizes de Coeficientes Tranpostas
awT = aw'; aeT = ae';
asT = as'; anT = an';
apT = ap';

% Termo fonte linearizado
SuT = reshape(Su', vc, 1);

% Matriz global do problema
M = zeros(vc);

for i=1:vc
    M(i,i) = apT(i);        % diagonal principal
    
    if i <= vc -1           % tri-diagonais
        M(i,i+1) = -aeT(i); 
        M(i+1,i) = -awT(i+1);
    end
    
    if i<= vc - nx            % penta-diagonais
        M(i,i + nx)= -asT(i);
        M(i + nx,i)= -anT(i+nx);
    end
end

% Solução do problema
T_inv = M^(-1)*SuT;

tempo_inv = toc;
%% TDMA
tic
% Coeficientes do TDMA
Aj = zeros(row,col);
Cj = zeros(row,col);
Clj= zeros(row,col);

% Incialização do campo da propriedade phi
phi = 100*zeros(row, col);

% Iniciando os residuos e iterações
ResG     = [];                  % Global
ResI     = zeros(row,col);      % Matriz de residuo local
iter     = 0; 
ResImax   = [];
 
% Tolerância do TDMA
tol = 1e-3;
CP  = 1; % Crit�rio de Parada 

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
%% Solução Analítica
T = zeros(row,col);

x_cont = col; %contador de VCs x
y_cont = row; %contador de VCs y

% Calcula a solucao analitica no centro de cada volume de controle
% a malha considerado na solucao numerica

for x = (dx/2:dx:L-dx/2) 
    for y = (dy/2:dy:W-dy/2)
        Sum = 0;
        for n = 1:100
            Sum = Sum + ( (-1)^(n+1) + 1 ) / n * sin(n*pi*x/L) * sinh(n*pi*y/L) / sinh(n*pi*W/L);
        end

        T(x_cont, y_cont) = 2*Sum/pi * (T2-T1) + T1;
        y_cont = y_cont-1;
    end
y_cont = row;
x_cont = x_cont -1;
end
%% Plot Solução pela inversa
Temp_fig = figure('Units','normalized');
imagesc(reshape(T_inv,row,col)')
colorbar;
axis ij;

title('Temperatura (�C) - Invers�o de matriz')
xlabel('x [m]')
ylabel('y [m]')

%% Plot Solução TDMA
Temp_fig = figure('Units','normalized');
imagesc(reshape(phi',row,col)')
colorbar;
axis ij;
title('Temperatura (�C) - TDMA');
xlabel('x [m]');
ylabel('y [m]');

Temp_fig = figure('Units','normalized');
plot(1:iter,ResG,'b')
title('Evolu��o do residuo Global');
xlabel('Itera��es');
ylabel('Residuo Global');
grid()

Temp_fig = figure('Units','normalized');
plot(1:iter,ResImax)
title('Evolu��o do m�ximo residuo local quadr�tico');
xlabel('Itera��es');
ylabel('Residuo Local ao quadrado');
grid()

%% Plot Solução Analítica
Temp_fig = figure('Units','normalized');
imagesc(reshape(T',row,col));
colorbar;
axis ij;
title('Solu��o anal�tica');
xlabel('x [m]');
ylabel('y [m]');
%% Plot Erro com Analítica
Temp_fig = figure('Units','normalized');
erro = reshape(T',row,col) - phi;
imagesc(erro)
colorbar;
axis ij;
title('Erro TDMA x Solu��o Anal�tica');

Temp_fig = figure('Units','normalized');
erro = reshape(T_inv,row,col)' - phi;
imagesc(erro)
colorbar;
axis ij;
title('Erro TDMA x Solu��o Inversa');








