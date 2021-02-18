clear
clc

// a)

// Método de la Potencia

// valor es el autovalor dominante
// zn es el autovector asociado

function [valor,zn] = m_potencia(A,z0,eps,maxIter)
    valor = 0;
    iter = 0;
    
    // iteración inicial
    w = A*z0;
    zn = w / norm(w,%inf); // podríamos usar cualquier norma
    [m,j] = max(abs(w)); // elijo la componente maxima de w y su índice
    valor = m / zn(j);
    iter = iter + 1;
    
    // iteraciones
    while (norm(z0-zn,%inf) > eps) && (iter < maxIter)
        z0 = zn;
        w = A*z0;
        zn = w / norm(w,%inf);
        [m,j] = max(abs(w));
        valor = m / zn(j);
        iter = iter + 1;
    end
    
    printf("Iteraciones: %d\n", iter);
endfunction

A1 = [6,4,4,1;4,6,1,4;4,1,6,4;1,4,4,6];
[lamb1,vect1] = m_potencia(A1,[1,0,0,0]',0.00001,1000);
printf("Matriz:\n");
disp(A1);
printf("\nVector:\n");
disp(vect1);
printf("\nAutovalor:\n");
disp(lamb1);
printf("\nVerificación:\tAutovalor:%f\n",max(abs(spec(A1))));
printf("A*vector:\n");
disp(A1*vect1);
printf("vector*lambda:\n");
disp(vect1*lamb1);

printf("---------------------------------------------------------------\n");

A2 = [12,1,3,4;1,-3,1,5;3,1,6,-2;4,5,-2,1];
[lamb2,vect2] = m_potencia(A2,[1,1,1,1]',0.00001,1000);
printf("Matriz:\n");
disp(A2);
printf("\nVector:\n");
disp(vect2);
printf("\nAutovalor:\n");
disp(lamb2);
printf("\nVerificación:\tAutovalor:%f\n",max(abs(spec(A2))));
printf("A*vector:\n");
disp(A2*vect2);
printf("vector*lambda:\n");
disp(vect2*lamb2);

printf("---------------------------------------------------------------\n");

//------------------------------------------------------------------------------//
function [valor,zn] = m_potenciab(A,z0,eps,maxIter)
    lambda = max(spec(A));
    
    valor = 0;
    iter = 0;
    
    // iteración inicial
    w = A*z0;
    zn = w / norm(w,%inf); // podríamos usar cualquier norma
    [m,j] = max(abs(w)); // elijo la componente maxima de w y su índice
    valor = m / zn(j);
    iter = iter + 1;
    err = abs(lambda-valor);
    
    // iteraciones
    while (err > %eps) // mientras que el error sea perceptible
        z0 = zn;
        w = A*z0;
        zn = w / norm(w,%inf);
        [m,j] = max(abs(w));
        valor = m / zn(j);
        iter = iter + 1;
        err = abs(lambda-valor);
    end
    
    printf("Iteraciones: %d\n", iter);
endfunction

[lamb2,vect2] = m_potenciab(A2,[1,1,1,1]',0.00001,1000);
printf("Matriz:\n");
disp(A2);
printf("\nVector:\n");
disp(vect2);
printf("\nAutovalor:\n");
disp(lamb2);
printf("\nVerificación:\tAutovalor:%f\n",max(abs(spec(A2))));
printf("A*vector:\n");
disp(A2*vect2);
printf("vector*lambda:\n");
disp(vect2*lamb2);
