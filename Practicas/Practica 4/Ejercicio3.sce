// a)
function [X,a] = gausselimsim(A,B)
    [nA,mA] = size(A);
    [nB,mB] = size(B);
    if nA <> mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif nA <> nB then
        error('gausselim - dimensiones incompatibles entre A y B');
        abort;
    end
    
    a = [A B]; // Matriz aumentada
    for i = 1:(nA-1) do // Eliminación progresiva
        for j = (i+1):nA do
           mji = a(j,i)/a(i,i); // multiplicador de fila
           a(j,i) = 0; // hacemos cero los elementos debajo de aii
           a(j,(i+1):(mA+mB)) = a(j,(i+1):(mA+mB)) - mji*a(i,(i+1):(mA+mB));
           // Cambiamos la fila j por la diferencia entre la fila j y la fila i multiplicada por mjk. mA+mB es la cantidad de columnas de la matriz aumentada [A B], y nA la cantidad de filas de la misma 
        end
    end
    // Sustitución regresiva
    X(nA,1:mB) = a(nA,(nA+1):(nA+mB))/a(nA,nA); // Calculamos la fila mB de X despejando de la ecuación mB-ésima de AX = B
    for i = (nA-1):-1:1 do
        X(i,1:mB) = (a(i,(mA+1):(mA+mB)) - (a(i,(i+1):mA)*X((i+1):mA,1:mB)))./a(i,i); // Calculamos la fila i de X despejando de la ecuación i-ésima de AX = B
    end
endfunction

// b)
A = [1,2,3;3,-2,1;4,2,-1];
B = [14,9,-2;2,-5,2;5,19,12];
[X,a] = gausselimsim(A,B);
disp(a);
disp(X);

// c)
// Queremos calcular la inversa de A, A-1, por eliminación gaussiana. Por definición de matriz inversa, AA-1 = A-1A = I, con I la matriz identidad. Luego, si en el algoritmo consideramos B = I, la solución de AX = B será X = A-1
printf("\nCálculo de la inversa:\n");
I = eye(3,3);
[X,a] = gausselimsim(A,I);
disp(a);
disp(X);
