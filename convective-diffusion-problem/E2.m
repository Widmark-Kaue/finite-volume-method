%% Head
clc; 
clear all; 
close all;
%format short eng;

%% Problem parameters

% Geometry
L  = 5e-3;           %[m]
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

%% Mesh

% Number of VC's
nx = 5; 
nr = 5; 
vc = nx*nr;      % total

% Length of VC's
dx = L/nx;
dr = R/nr;

row = nr; col = nx;

%% Coefficients matrix - UDS

%hvsd = @(x) [0*(x == 0) + (x > 0)]; %heaviside function (deg)
%
% rn : radial position of the N surfaces 
% rs : radial position of the S surfaces
% rm : radial position in vc's center.
% rm = [ rm_{1},...,rm_{i-1}, rm_{i}, rm_{i+1},..., rm_{n}]
% vc :                  N  ,   P   ,   S
%

rm  = R-dr/2:-dr :dr/2;
rn  = R     :-dr :dr;
rs  = R-dr  :-dr :0;

F   = rho.*u(rm);
D   = k/cp/dr;
G   = pi*(rn.^2 - rs.^2);

% Coefficients
ae  = zeros(row, col);

aw  = repmat(F', 1, col).*repmat(G',1,col);
aw (:,1) = 0;

as  = repmat(2*pi*D*dx , row, col).*repmat(rs',1,col);
as(row,:)= 0;

an  = repmat(D*2*pi*dx, row, col).*repmat(rn',1,col);
an(1,:)  = 0;

ap  = aw + an + as + ae;

%%  Boundary Conditions Aplication
Su  = zeros(row, col);

Su(:,1) = rho*um*T0*G;                  % west face 
ap(:,1) = ap(:,1) + (rho*um*G)';

Su(1,:) = Su(1,:) + k*R*q*2*pi*dx/cp;   % north face
                                        % Sp = 0
% south face: 
% Su = 0 ; Sp = 0

%% Solve - Matrix inversion method
M   = matrix_inv(ap,an,as,aw,ae);
SuL = reshape(Su', vc, 1);
T_inv = M^(-1)*SuL;


close all

T_2 = reshape(T_inv, col,row)'; 
T_2(row + 1:2*row,:) = flipud(T_2);
T_2 = T_2 - 273;

imagesc(T_2)
colorbar
%% Solve - TDMA method
[T_tdma, iter, Res, tempo] = tdma2d(ap,an,as,aw,ae,Su); 

close all

T_3 = T_tdma; 
T_3(row + 1:2*row,:) = flipud(T_3);
T_3 = T_3 - 273;

Tm_td = sum(T_tdma(:,end))/nr;
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
Tr = @(r) Ts - aux1*aux2(r) - 273;

%% Results

%Error
er_T_m  = Tm - Tm_td;
er_Ts   = Ts - T_tdma(1,end);

% Temperature profile
figure
hold on
r   = linspace(-R,R);
rm2 = R-dr/2: -dr : -(R-dr/2);
plot(Tr(r), r, 'b', 'LineWidth',2)
plot(T_3(:,end), rm2, 'ks')
legend('Analítica', 'Numérico')
grid on;
title('Perfil de Temperatura na Saída')
xlabel('Temperatura [K]')
ylabel('Posição Radial [m]' )
hold off

%%
% Temperature distribution
figure
r = -R:dr:R;
x = 0:dx:L;
subplot(2,1,1)
imagesc(x,r,T_3);
title('Distribuição de Temperaturas (TDMA) ')
xlabel('Quantidade Volumes Finitos \Deltax')
ylabel('Seção transversal do tubo [m]')
colorbar;

%{
subplot(2,1,2);
imagesc(T_3);
title('Distribuição de Temperaturas por TDMA')
xlabel('Quantidade Volumes Finitos \Deltax')
ylabel('Seção transversal do tubo [m]')
colorbar;
%}