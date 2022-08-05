%% Head
clc; 
clear all; 
close all;
format short eng;

%% Problem parameters

% Geometry
L  = 5e-3;          %[m]
Dd = 1e-3;           %[m]
R  = Dd/2;

% Fluid properties
k   = 1;            %[W/m.K]
rho = 100;          %[kg/m^2] 
mu  = 0.1;          %[Pa.s]
cp  = 1000;         %[J/kg.K]

% Flow properties 
Pr = 100;
um = 1;             %[m/s]

u  = @(r) 2*um*(1 - (r/R).^2); %perfil de velocidade

% Boundary conditions
T0 = 30 + 273;      %[K]     - temperatura de entrada do esc.
q  = 1000;          %[W/m^2] - fluxo de calor nas superficies

%% -------------------------------MFV-------------------------------------
%% Mesh

% Number of VC's
nx = [5 10 20 40]; 
nr = [5 10 20 40]; 
vc = nx.*nr;      % total

% Length of VC's
dx = L./nx;
dr = R./nr;

%% Coefficients matrix - UDS
%
% rn : radial position of the N surfaces 
% rs : radial position of the S surfaces
% rm : radial position in vc's center.
% rm = [ rm_{1},...,rm_{i-1}, rm_{i}, rm_{i+1},..., rm_{n}]
% vc :                  N  ,   P   ,   S
%
% south face: 
% Su = 0 ; Sp = 0

for i=1:length(nx)

    row = nr(i); col = nx(i);

    rm  = R-dr(i)/2 :-dr(i) :dr(i)/2;
    rn  = R         :-dr(i) :dr(i);
    rs  = R-dr(i)   :-dr(i) :0;

    F   = rho.*u(rm);           % convective flow
    D   = k/cp/dr(i);           % difusive flow

    G   = pi*(rn.^2 - rs.^2);
    Lr  = 2*pi*D*dx(i);


    % Coefficients
    ae{i}  = zeros(row, col);

    aw{i}  = repmat(F', 1, col).*repmat(G',1,col);
    aw{i}(:,1) = 0;

    as{i}  = repmat(Lr, row, col).*repmat(rs',1,col);
    as{i}(row,:)= 0;

    an{i} = repmat(Lr, row, col).*repmat(rn',1,col);
    an{i}(1,:)  = 0;

    ap{i} = aw{i} + an{i} + as{i} + ae{i};

    % Boundary Conditions Aplication
    Su{i} = zeros(row, col);

    Su{i}(:,1) = rho*um*T0*G;                           % west face 
    ap{i}(:,1) = ap{i}(:,1) + (rho*um*G)';

    Su{i}(1,:) = Su{i}(1,:) + k*R*q*2*pi*dx(i)/cp;      % north face
                                                        % Sp = 0                                                                                             
end
%% -------------------------------SOLVE------------------------------------
%% Solve - Matrix inversion method
for i=1:length(nx)
    row = nr(i); col = nx(i);
    tic
    M   = matrix_inv(ap{i},an{i},as{i},aw{i},ae{i});
    SuL = reshape(Su{i}', vc(i), 1);
    T_inv = M^(-1)*SuL;

    T_2{i} = reshape(T_inv, col,row)'; 
    T_2{i}(row + 1:2*row,:) = flipud(T_2{i});
    T_2{i} = T_2{i} - 273;
    t_inv(i) = toc;
end
%% Solve - TDMA method
for i=1:length(nx)
    row = nr(i); col = nx(i);
    [T_tdma, iter(i), Res{i}, tempo(i)] = tdma2d(ap{i},an{i},as{i},aw{i},ae{i},Su{i}); 

    T_num{i} = T_tdma; 
    T_num{i}(row + 1:2*row,:) = flipud(T_num{i});
    T_num{i} = T_num{i} - 273;

    Tm_td(i) = sum(T_tdma(:,end))/row;
    Ts_td(i) = T_tdma(1,end);
end
%% Solve - Analitic
m_dot = rho*um*pi*Dd^2/4;
P     = pi*Dd;
alpha = k/cp/rho;
dTm   = q*P/(m_dot*cp);

% Outlet
Tm = T0 + q*P/(m_dot*cp)*L;
Ts = Tm + (11/48)*q*Dd/k;

aux1 = 2*um*R^2*dTm/alpha;
aux2 = @(r) 3/16 + 1/16 *(r/R).^4 - 1/4 * (r/R).^2;
Trf  = @(r) Ts - aux1*aux2(r) - 273;

r  = linspace(-R,R);
Tr = Trf(r);

%% ------------------------------RESULTS----------------------------------
%% Error
er_Tm  = (Tm - Tm_td);
er_Ts  = (Ts - Ts_td);

figure
hold on
plot(vc, er_Tm,'b-o')
plot(vc, er_Ts,'r-o')
xlabel('Nº de Volumes de Controle')
ylabel('Erro')
legend('Temp. Média','Temp. da Superfície')
grid on
hold off

%% Temperature profile

figure
rm2 = R-dr(end)/2: -dr(end) : -(R-dr(end)/2);
hold on
plot(Tr, r.*1e3, 'b', 'LineWidth', 2)
plot(T_num{end}(:,end), rm2.*1e3, 'k--.')
legend('Analítica', 'Numérico')
hold off
grid on;
title(['(',num2str(nx(end)),'x',num2str(nx(end)),')'])
xlabel('Temperatura [ºC]')
ylabel('Posição Radial [mm]' )

%% Temperature distribution

hvsd = @(x) [0*(x == 0) + (x > 0)]; %heaviside function (deg)

for i=1:length(nx)
    j = mod(i,2)+ 2*mod(i-1,2);
    if mod(i,2)~=0
        figure
    end
    r = R:-dr(i):-R;
    x = 0:dx(i):L;
    subplot(2,1,j)
    imagesc(x.*1e3,r.*1e3,T_num{i});
    title(['Malha ',num2str(nx(i)),'x',num2str(nr(i))]);
    xlabel('Comprimento [mm]')
    ylabel('Seção transversal do tubo [mm]')
    colorbar;
end
%% Residues
for i=1:length(nx)
    figure
    hold on
    I = 1:iter(i);
    plot(I, Res{i}(1,:), 'k');
    plot(I, Res{i}(2,:), 'b');
    legend('Resíduo Local Quadrático Máximo', 'Residuo Normalizado')
    title(['Evolução dos Resíduos (',num2str(nx(i)),'x', num2str(nr(i)),')'])
    xlabel('Iterações')
    ylabel('Resíduos')
    grid on
    hold off
end
%% Error



