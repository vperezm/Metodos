clear
clc

// Método de Jacobi
// xi^(k+1) = 1/aii (bi - sum_(j=1,j!=i)^(n) aij xj^(k))

function [x,c] = jacobi(A,b,x0,eps)
    n = size(A,1);
    x = x0;
    xk = x;
    c = 0;
    
    // Iteración inicial
    for i = 1:n
        suma = 0;
        for j = 1:i-1
            suma = suma + A(i,j)*xk(j);
        end
        for j = i+1:n
            suma = suma + A(i,j)*xk(j);
        end
        x(i) = 1/A(i,i)*(b(i) - suma);
    end
    c = c + 1;
    
    // Iteraciones
    while (norm(x-xk) > eps) // no hace falta el valor absoluto porque norm(x) >= 0 forall x
        xk = x;
        for i = 1:n
            suma = 0;
            for j = 1:i-1
                suma = suma + A(i,j)*xk(j);
            end
            for j = i+1:n
                suma = suma + A(i,j)*xk(j);
            end
            x(i) = 1/A(i,i)*(b(i) - suma);
        end
        c = c + 1;
    end
endfunction

//------------------------------------------------------------------------------//
// Método de Gauss-Seidel
// xi^(k+1) = 1/aii (bi - sum_(j=1)^(i-1) aij xj^(k+1) - sum_(j=i+1)^(n) aij xj^(k))

function [x,c] = gauss_seidel(A,b,x0,eps)
    n = size(A,1);
    x = x0;
    xk = x;
    c = 0;
    
    // Iteración inicial
    for i = 1:n
        suma = 0;
        for j = 1:i-1
            suma = suma + A(i,j)*x(j);
        end
        for j = i+1:n
            suma = suma + A(i,j)*x(j);
        end
        x(i) = 1/A(i,i)*(b(i) - suma);
    end
    c = c + 1;
    
    // Iteraciones
    while (norm(x-xk) > eps) // no hace falta el valor absoluto porque norm(x) >= 0 forall x
        xk = x;
        for i = 1:n
            suma = 0;
            for j = 1:i-1
                suma = suma + A(i,j)*x(j);
            end
            for j = i+1:n
                suma = suma + A(i,j)*x(j);
            end
            x(i) = 1/A(i,i)*(b(i) - suma);
        end
        c = c + 1;
    end
endfunction
//------------------------------------------------------------------------------//

x0 = zeros(5,1);
eps = 10^(-6);

A = [10,1,2,3,4;1,9,-1,2,-3;2,-1,7,3,-5;3,2,3,12,-1;4,-3,-5,-1,15];
b = [12;-27;14;-17;12];

[x,c] = jacobi(A,b,x0,eps);
printf("Método de Jacobi luego de %d iteraciones:\n", c);
disp(x);
printf("\nVerificación:\n");
disp(A*x);
disp(b);
