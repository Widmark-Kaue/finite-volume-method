%% Head
clc; 
clear all; 
close all;
format short eng;

%% Problem parameters
% Geometry
W   = 0.5;              %[m]
L   = 0.5;              %[m]

% Fluid properties
alpha  = 1;             %[m^2/s]

% Boundary conditions
T1  = 0  ;      %[K]
T2  = 100;      %[K]

%% -------------------------------MFV-------------------------------------
%% Mesh
% Number of VC's
nx = 10; 
ny = 10; 
vc = nx*ny; % total

% Length of VC's
dx = L/nx;
dy = W/ny;

%% Time Step
endTime = 10;          %[s] 
dt      = 0.1;           %[s]

theta   = 1;
if theta ~=1
    aux1 = dx*dy;
    aux2 = (dy/dx + dx/dy);
    aux3 = (1-theta);
    aux4 = 2*alpha;
    dt = 0.5*aux1/aux2/aux3/aux4;
end

%% Coefficients matrix - CDS

row = ny; col = nx;

aw = repmat(dy/dx, row, col); 
aw(:,1)     = 0;

ae = repmat(dy/dx, row, col);
ae(:,col)   = 0;

as = repmat(dx/dy, row, col);
as(row,:)   = 0;

an = repmat(dx/dy, row, col);
an(1,:)     = 0;

ap0 = repmat(dx*dy/alpha/dt, row, col);
ap0(1,:)    = ap0(1,:)/2;
ap0(row,:)  = ap0(row,:)/2;
ap0(:,1)    = ap0(:,1)/2;
ap0(:,col)  = ap0(:,col)/2;

ap = (aw + ae + as + an)*theta + ap0;
Ap0= ap0 - (aw + ae + as + an).*(1-theta);

%Boundary Conditions Aplication
Su = zeros(row, col);

Su(1,:)     = 2*T2*dx/dy;
Su(row,:)   = Su(row,:)+ 2*T1*dx/dy;
Su(:,1)     = Su(:,1)  + 2*T1*dy/dx;
Su(:,col)   = Su(:,col)+ 2*T1*dy/dx;
                                
                       %-Sp
ap(1,:)     = ap(1,:)  + 2*dx/dy * theta;
ap(row,:)   = ap(row,:)+ 2*dx/dy * theta;
ap(:,1)     = ap(:,1)  + 2*dy/dx * theta;
ap(:,col)   = ap(:,col)+ 2*dy/dx * theta;

%Packed
A   = struct('ap',ap,'an',an,'as',as,'aw',aw, 'ae', ae);
Tt  = struct('Ap0',Ap0,'theta',theta); 

%% -------------------------------SOLVE------------------------------------
%% Solve - TDMA method
T_td = [];
for t=1:endTime/dt
    t
    if t == 1
        T_td(:,:,t) = tdma2d(A,Su);
    else
        T_td(:,:,t) = tdma2d(A,Su,T_td(:,:,t-1),Tt);
        
        res = max(max((T_td(:,:,t) - T_td(:,:,t-1)).^2));
        %{
            if res < 1e-3
            break
            end
        %}
    end
end
%% Solve - Analitic
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
%% ------------------------------RESULTS----------------------------------
%% Plot Solução pela inversa
Temp_fig = figure('Units','normalized');
imagesc(reshape(T_inv,row,col)')
colorbar;
axis ij;

title('Temperatura (ºC) - Inversão de matriz')
xlabel('x [m]')
ylabel('y [m]')

%% Teste plot
tam = size(T_td);
T2  = T_td;
x   = 0:nx/2:W*1e2;
y   = 0:ny/2:L*1e2; 
figure
for i=1:tam(end)
    imagesc(x,y,T2(:,:,i));
    title(['Distribuição de temperatura (t = ', num2str(i), 's )'])
    xlabel('Largura [cm]')
    ylabel('Altura [cm]')
    colorbar;
    anima(i) = getframe;
end
%%
%% Plot Solução TDMA
Temp_fig = figure('Units','normalized');
imagesc(reshape(phi',row,col)')
colorbar;
axis ij;
title('Temperatura (ºC) - TDMA');
xlabel('x [m]');
ylabel('y [m]');

Temp_fig = figure('Units','normalized');
plot(1:iter,ResG,'b')
title('Evolução do residuo Global');
xlabel('Iterações');
ylabel('Residuo Global');
grid()

Temp_fig = figure('Units','normalized');
plot(1:iter,ResImax)
title('Evolução do máximo residuo local quadrático');
xlabel('Iterações');
ylabel('Residuo Local ao quadrado');
grid()

%% Plot Solução Analítica
Temp_fig = figure('Units','normalized');
imagesc(reshape(T',row,col));
colorbar;
axis ij;
title('Solução analítica');
xlabel('x [m]');
ylabel('y [m]');
%% Plot Erro com Analítica
Temp_fig = figure('Units','normalized');
erro = reshape(T',row,col) - phi;
imagesc(erro)
colorbar;
axis ij;
title('Erro TDMA x Solução Analítica');

Temp_fig = figure('Units','normalized');
erro = reshape(T_inv,row,col)' - phi;
imagesc(erro)
colorbar;
axis ij;
title('Erro TDMA x Solução Inversa');








