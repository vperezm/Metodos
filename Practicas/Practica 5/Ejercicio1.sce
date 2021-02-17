clear
// a) y b)

// Matrices N para el método de Jacobi y Gauss-Seidel
function [Nj,Ng] = matriz_N(A)
    n = size(A,1);
    // inicializamos las matrices en cero
    Nj = zeros(n,n);
    Ng = zeros(n,n);
    // completamos
    for i = 1:n
        Nj(i,i) = A(i,i);
        for j = 1:i
            Ng(i,j) = A(i,j);
        end
    end
endfunction

I = eye(3,3);

//------------------------------------------------------------------------------//
printf("---------------------------------------------------------------------------");
A1 = [0,2,4;1,-1,-1;1,-1,2];
b1 = [0;0.375;0];
// A1 no es diagonal dominante  
[Nj1,Ng1] = matriz_N(A1);
disp(Nj1);
printf("det(Nj) = %f\n", det(Nj1)); // = 0 -> Nj1 no es inversible
disp(Ng1);
printf("det(Ng) = %f\n\n", det(Ng1)); // = 0 -> Ng1 no es inversible

printf("det(A) = %f\n", det(A1)); // = -6 -> A1 es no singular
// puedo hacer un intercambio de filas para lograr Nj1 y Ng1 no singulares

// A(1) <-> A(2)
// b(1) <-> b(2)
printf("\nCambio de filas: A(1) <-> A(2), b(1) <-> b(2)\n");
A1 = [1,-1,-1;0,2,4;1,-1,2];
b1 = [0.375;0;0];  
[Nj1,Ng1] = matriz_N(A1);
disp(Nj1);
printf("det(Nj) = %f\n", det(Nj1)); // = 4 -> Nj1 es inversible
disp(Ng1);
printf("det(Ng) = %f\n\n", det(Ng1)); // = 4 -> Ng1 es inversible

// Veo si puedo usar Teorema 1 para asegurar la convergencia
disp(norm(I - inv(Nj1)*A1)); // >= 1 -> no es suficiente
disp(norm(I - inv(Ng1)*A1)); // >= 1 -> no es suficiente

// Veo si puedo usar Corolario 1 para asegurar la convergencia
disp(max(abs(spec(I - inv(Nj1)*A1)))); // >= 1 -> no es suficiente 
disp(max(abs(spec(I - inv(Ng1)*A1)))); // >= 1 -> no es suficiente

//------------------------------------------------------------------------------//
printf("---------------------------------------------------------------------------");
A2 = [1,-1,0;-1,2,-1;0,-1,1.1];
b2 = [0;1;0];
// A2 no es diagonal dominante
[Nj2,Ng2] = matriz_N(A2);
disp(Nj2);
printf("det(Nj) = %f\n", det(Nj2)); // = 2 -> Nj2 es inversible
disp(Ng2);
printf("det(Ng) = %f\n\n", det(Ng2)); // = 2 -> Ng2 es inversible

// Veo si puedo usar Teorema 1 para asegurar la convergencia
disp(norm(I - inv(Nj2)*A2)); // >= 1 -> no es suficiente
disp(norm(I - inv(Ng2)*A2)); // >= 1 -> no es suficiente

// Veo si puedo usar Corolario 1 para asegurar la convergencia
disp(max(abs(spec(I - inv(Nj2)*A2)))); // = 0.9770084 < 1 -> converge
disp(max(abs(spec(I - inv(Ng2)*A2)))); // = 0.9545455 < 1 -> converge
