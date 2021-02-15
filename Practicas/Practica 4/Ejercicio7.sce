// Factorización PA = LU a partir de la eliminación de Gauss con pivoteo parcial

function [U,L,P] = palu_gausspp(A)
    [n,m] = size(A);
    if n <> m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    
    U = A;
    L = eye(n,n);
    P = eye(n,n);
    
    for k = 1:n-1
        // Pivoteo
        // kpivot = i que maximiza |Uik|
        kpivot = k; amax = abs(U(k,k));
        for i=k+1:n
            if abs(U(i,k))>amax then
                kpivot = i; amax = U(i,k);
            end
        end
        
        // U(k,k:n) <-> U(i,k:n)
        temp = U(k,k:n);
        U(k,k:n) = U(i,k:n);
        U(i,k:n) = temp;
        
        // L(k,1:k-1) <-> L(i,1:k-1)
        temp = L(k,1:k-1);
        L(k,1:k-1) = L(i,1:k-1);
        L(i,1:k-1) = temp;
        
        // P(k) <-> P(i)
        temp = P(k);
        P(k) = P(i);
        P(i) = temp;
        
        for j = k+1:n
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
        end
    end
endfunction

A = [2,1,1,0;4,3,3,1;8,7,9,5;6,7,9,8];
[U,L,P] = palu_gausspp(A);
disp(U);
disp(L);
disp(P);
