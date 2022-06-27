%% Solução pela inversa
function M = matrix_inv(ap, an, as, aw, ae)

[nr, nx] = size(ap);

% Matrizes de Coeficientes Tranpostas
awT = aw'; aeT = ae';
asT = as'; anT = an';
apT = ap';

% Matriz global do problema

M = diag(apT(1:end));
M = M - diag(aeT(1:end-1), 1) - diag(awT(2:end), -1);
M = M - diag(asT(1:end-nx), nx) - diag(anT(nx+1:end), -nx);

%{
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
%}
end