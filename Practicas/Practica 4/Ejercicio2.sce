//----------------------------------------------------------------------------------//
/*
// c)
// Función gausselim modificada para contar el numero de operaciones realizadas
function [x,a,c] = gausselim(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

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

c = 0; // Variable para contar las operaciones realizadas

// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            c = c + 3;
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n);
c = c+1;
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k);
        c = c + 2;
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
    c = c + 2;
end;
endfunction
*/

//----------------------------------------------------------------------------------//
// d)
// Función gausselim modificada para usar submatrices en Scilab
function [x,a] = gausselim(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

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

// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        a(i,k+1:n+1) = a(i,k+1:n+1) - a(k,k+1:n+1)*a(i,k)/a(k,k);
        a(i,1:k) = 0;  // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    x(i) = (a(i,n+1) - a(i,i+1:n)*x(i+1:n))/a(i,i);
end;
endfunction
//----------------------------------------------------------------------------------//
// b)

// i)
printf("-------------------------------------\n");
A1 = [1,1,0,3;2,1,-1,1;3,-1,-1,2;-1,2,3,-1];
b1 = [4;1;-3;4];
[x1,a1] = gausselim(A1,b1);
disp(x1);
disp(a1);

// ii)
printf("-------------------------------------\n");
A2 = [1,-1,2,-1;2,-2,3,-3;1,1,1,0;1,-1,4,3];
b2 = [-8;-20;-2;4];
[x2,a2] = gausselim(A2,b2);
disp(x2);
disp(a2);

// iii)
printf("-------------------------------------\n");
A3 = [1,1,0,4;2,1,-1,1;4,-1,-2,2;3,-1,-1,2];
b3 = [2;1;0;-3];
[x3,a3] = gausselim(A3,b3);
disp(x3);
disp(a3);
