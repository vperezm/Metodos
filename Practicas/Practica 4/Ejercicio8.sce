// Factorización PA = LU a partir de la eliminación de Gauss con pivoteo parcial

function [L,U,P] = palu_gausspp(A)
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

// a)
A1 = [1.012,-2.132,3.104;-2.132,4.096,-7.013;3.104,-7.013,0.014];

[L1,U1,P] = palu_gausspp(A1);
disp(L1);
disp(U1);

[L1s,U1s] = lu(A1);
disp(L1s);
disp(U1s);

// b)
A2 = [2.1756,4.0231,-2.1732,5.1967;-4.0231,6.0000,0,1.1973;-1.0000,5.2107,1.1111,0;6.0235,7.0000,0,4.1561];
[L2,U2,P] = palu_gausspp(A2);
disp(L2);
disp(U2);

[L2s,U2s] = lu(A2);
disp(L2s);
disp(U2s);
