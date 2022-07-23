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
nx = 20; 
ny = 20; 
vc = nx*ny; % total

% Length of VC's
dx = L/nx;
dy = W/ny;

%% Time Step
endTime = 0.1;          %[s] 
dt0     = 0.001;        %[s]

theta   = [1 0.5 0]; 
% theta = 
%   1     - fully impicit scheme
%   0.5   - Crank-Nicolson scheme
%   0     - Explicit scheme 


% dt for theta != 1
aux1 = dx*dy;
aux2 = (dy/dx + dx/dy);
aux3 = @(th) (1-th);
aux4 = 2*alpha;

fdt = @(th) 0.5*aux1/aux2/aux3(th)/aux4;


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


ap0 = {};
ap  = {};
Ap0 = {};
for i=1:length(theta)
    
    if theta(i)~=1
        ap0{i} = repmat(dx*dy/alpha/fdt(theta(i)), row, col);
    else
        ap0{i} = repmat(dx*dy/alpha/dt0, row, col);
    end
    ap0{i}(1,:)    = ap0{i}(1,:)/2;
    ap0{i}(row,:)  = ap0{i}(row,:)/2;
    ap0{i}(:,1)    = ap0{i}(:,1)/2;
    ap0{i}(:,col)  = ap0{i}(:,col)/2;
    
    
    ap{i} = (aw + ae + as + an)*theta(i) + ap0{i};
    Ap0{i}= ap0{i} - (aw + ae + as + an).*(1-theta(i));
    
                                 %-Sp
    ap{i}(1,:)     = ap{i}(1,:)  + 2*dx/dy * theta(i);
    ap{i}(row,:)   = ap{i}(row,:)+ 2*dx/dy * theta(i);
    ap{i}(:,1)     = ap{i}(:,1)  + 2*dy/dx * theta(i);
    ap{i}(:,col)   = ap{i}(:,col)+ 2*dy/dx * theta(i);
end

%Boundary Conditions Aplication
Su = zeros(row, col);

Su(1,:)     = 2*T2*dx/dy;
Su(row,:)   = Su(row,:)+ 2*T1*dx/dy;
Su(:,1)     = Su(:,1)  + 2*T1*dy/dx;
Su(:,col)   = Su(:,col)+ 2*T1*dy/dx;
               

%% -------------------------------SOLVE------------------------------------
%% Solve - TDMA method
T_td        = {};
tol_tr      = 1e-3;
tempo       = [];

for i=1:length(theta)
    tic
    %Packed
    A   = struct('ap',ap{i},'an',an,'as',as,'aw',aw, 'ae', ae);
    Tt  = struct('Ap0',Ap0{i},'theta',theta(i)); 
   
    % TDMA iterative initial
    fprintf('################### THETA = %f ####################\n', theta(i))
    if theta(i) ~= 1
        dt = fdt(theta(i));
    else
        dt = dt0;
    end
    sq_res      = [1];
    res         = [1]; 
    for t=1:endTime/dt
        fprintf('---------------- Step: %d ----------------\n', t)
        fprintf('Step time: %f\n',dt*t);
        if t == 1
            [T_td{i}(:,:,t),iter] = tdma2d(A,Su);
        else
            [T_td{i}(:,:,t),iter] = tdma2d(A,Su,T_td{i}(:,:,t-1),Tt);

            sq_res(t) = max(max((T_td{i}(:,:,t) - T_td{i}(:,:,t-1)).^2));
            res(t)    = sum(sum(abs(T_td{i}(:,:,t) - T_td{i}(:,:,t-1))));
            if sq_res(t) < tol_tr
                fprintf('CONVERGIU!')
                break
            end
        end
        fprintf('Iterarions: %d\n',iter);
        fprintf('Residue: %f\n', sq_res(t));
    end
    tempo(i) = toc;
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
%% plot temperature field - theta = 1
tam = size(T_td{1});
T2  = T_td{1};
x   = 0:nx/2:W*1e2;
y   = 0:ny/2:L*1e2; 
for i=1:tam(end)
    figure
    imagesc(x,y,T2(:,:,i));
    title(['Distribuição de temperatura (t = ', num2str(i*dt0), 's ) - \theta = 1'])
    xlabel('Largura [cm]')
    ylabel('Altura [cm]')
    colorbar;
    anima(i) = getframe;
end

%% plot temperature field - theta = 0.5
tam = size(T_td{2});
T2  = T_td{2};
x   = 0:nx/2:W*1e2;
y   = 0:ny/2:L*1e2; 
figure
for i=1:tam(end)
    imagesc(x,y,T2(:,:,i));
    title(['Distribuição de temperatura (t = ', num2str(i*dt0), 's ) - \theta = 0.5'])
    xlabel('Largura [cm]')
    ylabel('Altura [cm]')
    colorbar;
    anima(i) = getframe;
end
%% plot temperature field - theta = 0
tam = size(T_td{3});
T2  = T_td{3};
x   = 0:nx/2:W*1e2;
y   = 0:ny/2:L*1e2; 
figure
for i=1:tam(end)
    imagesc(x,y,T2(:,:,i));
    title(['Distribuição de temperatura (t = ', num2str(i*dt0), 's ) - \theta = 0'])
    xlabel('Largura [cm]')
    ylabel('Altura [cm]')
    colorbar;
    anima(i) = getframe;
end
%% plot convergence

forms = ['k', 'b', 'r']
figure
hold on
for i=1:length(theta)
    tam = size(T_td{i});
    T3  = T_td{i};
    
    x = [];
    y = [];
    if i==1
        dt = dt0;
    else
        dt = fdt(theta(i));
    end
    for j=1:tam(end)
        y(j) = sum(T3(1,:,j))/nx;
        x(j) = j*dt;
    end
    plot(x, y, forms(i), 'LineWidth',2)
end
legend('\theta = 1', '\theta = 0.5', '\theta = 0')
title('Temperatura média na superfície superior')
xlabel('Tempo [s]')
ylabel('Temperatura [ºC]')
grid
hold off

%% plot convergence -  theta = 1
figure

tam = size(T_td{1});
T3  = T_td{1};

a = sum(T(1,:));
A = a*ones(1,length(tam));
x = [];
y = [];
for j=1:tam(end)
    y(j) = sum(T3(1,:,j))/nx;
    x(j) = j*dt;
end

plot(x,A,'r--','LineWidth',2)
plot(x, y, 'k', 'LineWidth',2)
legend('\theta = 1')
title('Temperatura média na superfície superior')
xlabel('Tempo [s]')
ylabel('Temperatura [ºC]')
grid
hold off
