clear
clc
// Agrego un máximo de iteraciones a los métodos, ya que por a) y b) se que
// no puedo asegurar la convergencia
// ESTO ESTA MAL. SI NO CONVERGE NO CONVERGE

//------------------------------------------------------------------------------//
// Método de Jacobi
// xi^(k+1) = 1/aii (bi - sum_(j=1,j!=i)^(n) aij xj^(k))

function x = jacobi(A,b,x0,eps,max_iter)
    n = size(A,1);
    x = x0;
    iter = 0;
    
    // Iteración inicial
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
    iter = iter + 1;
    
    // Iteraciones
    while (norm(x-xk) > eps) && (iter < max_iter)
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
        iter = iter + 1;
    end
    if iter == max_iter then
        error("No se pudo llegar a un resultado con el máximo de iteraciones");
        abort;
    end
endfunction

//------------------------------------------------------------------------------//
// Método de Gauss-Seidel
// xi^(k+1) = 1/aii (bi - sum_(j=1)^(i-1) aij xj^(k+1) - sum_(j=i+1)^(n) aij xj^(k))

function x = gauss_seidel(A,b,x0,eps,max_iter)
    n = size(A,1);
    x = x0;
    iter = 0;
    
    // Iteración inicial
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
    iter = iter + 1;
    
    // Iteraciones
    while (norm(x-xk) > eps) && (iter < max_iter)
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
        iter = iter + 1;
    end
    if iter == max_iter then
        error("No se pudo llegar a un resultado con el máximo de iteraciones");
        abort;
    end
endfunction

//------------------------------------------------------------------------------//
x0 = zeros(3,1);
eps = 10^(-2);
max_iter = 100000;

A1 = [1,-1,-1;0,2,4;1,-1,2];
b1 = [0.375;0;0];  
//printf("Método de Jacobi:\n");
//xj1 = jacobi(A1,b1,x0,eps,max_iter);
//disp(xj1);
//printf("Método de Gauss-Seidel:\n");
//xg1 = gauss_seidel(A1,b1,x0,eps,max_iter);
//disp(xg1);

// Ambos métodos fallaron para el sistema 1


printf("---------------------------------------------------------------------------\n");


A2 = [1,-1,0;-1,2,-1;0,-1,1.1];
b2 = [0;1;0];
printf("Método de Jacobi:\n");
xj2 = jacobi(A2,b2,x0,eps,max_iter);
disp(xj2);
printf("Verificación:\n");
disp(A2*xj2);
disp(b2);
printf("\nMétodo de Gauss-Seidel:\n");
xg2 = gauss_seidel(A2,b2,x0,eps,max_iter);
disp(xg2);
printf("Verificación:\n");
disp(A2*xg2);
disp(b2);
