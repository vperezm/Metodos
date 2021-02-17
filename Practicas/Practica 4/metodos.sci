//-----------------------------------------------------------------------------//
// SUSTITUCION REGRESIVA Y PROGRESIVA

// Sistemas triangulares superiores
function x = sust_reg(A,b)
    n = size(A,1);
    x(n) = b(n)/A(n,n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i,i+1:n)*x(i+1:n))/A(i,i);
    end
endfunction

// Sistemas triangulares inferiores
function x = sust_prog(A,b)
    n = size(A,1);
    x(1) = b(1)/A(1,1);
    for i = 2:n
        x(i) = (b(i) - A(i,1:i-1)*x(1:i-1))/A(i,i);
    end
endfunction

//-----------------------------------------------------------------------------//
// ELIMINACIÓN DE GAUSS

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

//-----------------------------------------------------------------------------//
// ELIMINACIÓN DE GAUSS PARA MÚLTIPLES SISTEMAS DE ECUACIONES

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

//-----------------------------------------------------------------------------//
// CÁLCULO DEL DETERMINANTE MEDIANTE EL MÉTODO DE GAUSS

function y = det_gauss(A)
    [nA,mA] = size(A);
    n = nA;
    a = A;
    for k=1:n
        for i=k+1:n
            a(i,k+1:n) = a(i,k+1:n) - a(k,k+1:n)*a(i,k)/a(k,k);
            a(i,1:k) = 0;  // no hace falta para calcular la solución x
        end;
    end;
    disp(a);
    y = 1;
    for i = 1:n
        y = y * a(i,i);
    end
endfunction

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
// MÉTODO DE GAUSS PARA MATRICES TRIDIAGONALES

// Con pivoteo parcial
function [x,a] = gausselimPPtri(A,b)
    [nA,mA] = size(A);
    [nb,mb] = size(b);
    if nA <> mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA <> nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end
    
    a = [A b]; // Matriz aumentada
    n = nA; // Tamaño de la matriz
    contador = 0;
    // Eliminación progresiva con pivoteo parcial
    for k = 1:n-1 do
        kpivot = k; amax = abs(a(k,k)); // pivoteo
        if abs(a(k+1,k)) > amax then
            temp = a(k+1,:);
            a(k+1,:) = a(k,:);
            a(k,:) = temp; // Intercambiamos las filas k y k+1
        end
        for j = k+1:n+1 do
            a(k+1,j) = a(k+1,j) - a(k,j)*a(k+1,k)/a(k,k);
            contador = contador + 3;
        end
    end
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    x(n-1) = (a(n-1,n+1)-a(n-1,n)*x(n))/a(n-1,n-1);
    contador = contador + 4;
    for i = n-2:-1:1 do
        x(i) = (a(i,n+1)-a(i,i+1)*x(i+1)-a(i,i+2)*x(i+2))/a(i,i);
        contador = contador + 5;
    end
    disp(contador);
endfunction

// Sin pivoteo
function [x,a] = gausselimtri(A,b)
    [nA,mA] = size(A);
    [nb,mb] = size(b);
    if nA <> mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA <> nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end
    
    a = [A b]; // Matriz aumentada
    n = nA; // Tamaño de la matriz
    contador = 0;
    // Eliminación progresiva sin pivoteo
    for k = 1:n-1 do
        mkk = a(k+1,k)/a(k,k);
        a(k+1,k+1) = a(k+1,k+1) - mkk*a(k,k+1);
        a(k+1,n+1) = a(k+1,n+1) - mkk*a(k,n+1);
        contador = contador + 5;
    end
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    contador = contador + 1;
    for i = n-1:-1:1 do
        x(i) = (a(i,n+1)-a(i,i+1)*x(i+1))/a(i,i);
        contador = contador + 3;
    end
    disp(contador);
endfunction

//-----------------------------------------------------------------------------//
// Factorización PA = LU a partir de la eliminación de Gauss con pivoteo parcial

function [U,L,P] = palu_gausspp(A)
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

//-----------------------------------------------------------------------------//
// MÉTODO DE DOOLITTLE

function [L,U] = doolittle(A)
    [n m] = size(A);
    L = eye(n,n);
    U = zeros(n,n);
    if n <> m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    for i = 1:n
        
        suma = 0;
        for k = 1:i-1
            suma = suma + L(i,k)*U(k,i:n);
        end
        U(i,i:n) = A(i,i:n) - suma;
        
        suma = 0;
        for k = 1:i-1
            suma = suma + L(i:n,k)*U(k,i);
        end
        L(i:n,i) = (A(i:n,i) - suma)/U(i,i);
    end
endfunction

//-----------------------------------------------------------------------------//
// FACTORIZACIÓN DE CHOLESKY

function [U, ind] = Cholesky(A)
eps = 1.0e-8
n = size(A,1)
U = zeros(n,n)
for k = 1:n
    if k==1 then
            t = A(k,k)
    else 
            t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    end

    if t <= eps
        printf("Matriz no definida positiva.\n")
        ind = 0
        return
    end
    U(k,k)= sqrt(t)
    for j = k+1:n
        if k==1 then 
                    U(k,j) = A(k,j)/U(k,k)
        else 
                    U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
        end
    end
end
ind = 1
endfunction

function x = chol_solver(A,b)
    [R,ind] = Cholesky(A);
    if ind == 0 then
        error('No se puede calcular');
        abort;
    end
    g = sust_prog(R',b);
    x = sust_reg(R,g);
endfunction
