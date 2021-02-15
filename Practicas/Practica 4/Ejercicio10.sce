clear
function [L,U] = doolittle(A)
    [n m] = size(A);
    L = eye(n,n);
    U = zeros(n,n);
    if n <> m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    for i = 1:n
        
        suma = 0;
        for k = 1:i-1
            suma = suma + L(i,k)*U(k,i:n);
        end
        U(i,i:n) = A(i,i:n) - suma;
        
        suma = 0;
        for k = 1:i-1
            suma = suma + L(i:n,k)*U(k,i);
        end
        L(i:n,i) = (A(i:n,i) - suma)/U(i,i);
    end
endfunction

A = [2,-1,-2;-4,6,3;-4,-2,8];
[L U] = doolittle(A);
disp(L);
disp(U);
