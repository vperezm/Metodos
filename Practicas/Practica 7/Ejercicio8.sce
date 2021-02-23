xdel
clear
clc

// Aproximación de mínimos cuadrados
// x,y vectores tales que (x(i),y(i)) son los puntos dados para realizar la aproximación, n es el grado del polinomio, P es el polinomio de aproximación de mínimos cuadrados, y err es el error cometido

function [P,err] = minimoscuadradosn(x,y,n)
    m = length(x);
    A = zeros(m,n+1); // A contiene los valores de las funciones Φ0(xi),...,Φn(xi) en la fila i-ésima, por lo que tiene n+1 columnas
    
    // Definimos la matriz A
    for i = 1:m
        for j = 1:n+1
            A(i,j) = x(i)^(j-1);
        end
    end
    
    // Calculamos AT A y AT b
    Amin = A'*A;
    bmin = A'*y';
    
    // Calculamos la solución 'a'
    a = inv(Amin)*bmin;
    
    // Definimos el polinomio de aproximación
    P = poly(a,"s",["coeff"]);
    
    // Calculamos el error de aproximación
    err = norm(A*a-y');
endfunction

x = [4,4.2,4.5,4.7,5.1,5.5,5.9,6.3,6.8,7.1];
y = [102.56,113.18,130.11,142.05,167.53,195.14,224.87,256.73,299.5,326.72];

printf("----------------------------------\n");
[p1,e1] = minimoscuadradosn(x,y,1);
printf("Grado 1:\n\n");
printf("Polinomio:\n");
disp(p1);
printf("Error:\n");
disp(e1);

printf("----------------------------------\n");
[p2,e2] = minimoscuadradosn(x,y,2);
printf("Grado 2:\n\n");
printf("Polinomio:\n");
disp(p2);
printf("Error:\n");
disp(e2);

printf("----------------------------------\n");
[p3,e3] = minimoscuadradosn(x,y,3);
printf("Grado 3:\n\n");
printf("Polinomio:\n");
disp(p3);
printf("Error:\n");
disp(e3);

// Gráfico

scatter(x,y,"fill");
//plot(x,horner(p1,x));
//plot(x,horner(p2,x));
plot(x,horner(p3,x));
