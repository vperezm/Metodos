clear
clc
//-----------------------------------------------------------------------------//
// MÉTODO DE GAUSS CON PIVOTEO PARCIAL

function [x,a] = gausselimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada
n = nA;    // Tamaño de la matriz

// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(a(k,k));  //pivoteo
    for i=k+1:n
        if abs(a(i,k))>amax then
            kpivot = i; amax = a(i,k);
        end;
    end;
    temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k);
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
end;
endfunction

//-----------------------------------------------------------------------------//
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

//-----------------------------------------------------------------------------//
// Función para crear la matriz pedida

function [A,b] = make_matrix(N)
    A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + ...
        diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
    b = ones(N,1);
endfunction

// Testing para la función make_matrix
/*
[A,b] = make_matrix(5);
disp(A);
disp(b);
*/

//-----------------------------------------------------------------------------//

N = 100;
printf("\nN = %d\n",N);
[A,b] = make_matrix(N);

tic();
[x,a] = gausselimPP(A,b);
t = toc();
printf("Método de eliminación de Gauss con pivoteo parcial\nt = %f\n", t);

x0 = zeros(N,1);
eps = 10^(-6);
tic();
[x,a] = gauss_seidel(A,b,x0,eps);
t = toc();
printf("Método de Gauss-Seidel con eps = 10^(-6)\nt = %f\n",t);

x0 = zeros(N,1);
eps = 10^(-11);
tic();
[x,a] = gauss_seidel(A,b,x0,eps);
t = toc();
printf("Método de Gauss-Seidel con eps = 10^(-11)\nt = %f\n",t);

//-----------------------------------------------------------------------------//

N = 500;
printf("\nN = %d\n",N);
[A,b] = make_matrix(N);

tic();
[x,a] = gausselimPP(A,b);
t = toc();
printf("Método de eliminación de Gauss con pivoteo parcial\nt = %f\n", t);

x0 = zeros(N,1);
eps = 10^(-6);
tic();
[x,a] = gauss_seidel(A,b,x0,eps);
t = toc();
printf("Método de Gauss-Seidel con eps = 10^(-6)\nt = %f\n",t);

x0 = zeros(N,1);
eps = 10^(-11);
tic();
[x,a] = gauss_seidel(A,b,x0,eps);
t = toc();
printf("Método de Gauss-Seidel con eps = 10^(-11)\nt = %f\n",t);
