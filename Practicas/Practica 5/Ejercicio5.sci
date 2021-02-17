clear
clc

// Calcula el w óptimo para A definida positiva y tridiagonal
function w = make_w(A)
    deff("re = p(T)","re = max(abs(spec(T)))");
    Tj = make_J(A);
    w = 2 / (1 + sqrt(1 - p(Tj)^2));
endfunction

// Calcula la matriz de Jacobi
function J = make_J(A)
    n = size(A,1);
    D = zeros(n,n);
    for i = 1:n
        D(i,i) = A(i,i);
    end
    I = eye(n,n);
    J = I - inv(D)*A;
endfunction

// Método de Sobrerrelajación (SOR)
function [x,c] = SOR(A,b,x0,eps,w)
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
        x(i) = (1-w)*x(i) + (w/A(i,i))*(b(i)-suma);
    end
    c = c + 1;
    
    // Iteraciones
    while (norm(x-xk) > eps)
        xk = x;
        for i = 1:n
            suma = 0;
            for j = 1:i-1
                suma = suma + A(i,j)*x(j);
            end
            for j = i+1:n
                suma = suma + A(i,j)*x(j);
            end
            x(i) = (1-w)*x(i) + (w/A(i,i))*(b(i)-suma);
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

A = [4,3,0;3,4,-1;0,-1,4];
b = [24;30;-24];

sol = [3;4;-5];

x0 = zeros(3,1);
eps = 10^(-7);

printf("\nMétodo de Gauss-Seidel:\n");
[x,c] = gauss_seidel(A,b,x0,eps);
disp(x);
disp(sol);
printf("Iteraciones: %d\n",c);

printf("\nMétodo SOR con parámetro de relajación óptimo:\n");
w = make_w(A);
[x,c] = SOR(A,b,x0,eps,w);
disp(x);
disp(sol);
printf("Iteraciones: %d\n",c);

