clear

// Funcion auxiliar para comprobar las dimensiones de las matrices ingresadas
function y = sz_check(A,b)
    szA = size(A);
    szb = size(b)
    if (szA(1) == szA(2)) && (szA(1) == szb(1)) && (szb(2) == 1) then
        y = szA(1);
    else
        y = 0;
    end
endfunction

// Función auxiliar para comprobar que la matriz ingresada es triangular superior
function y = tsup_check(A)
    y = 1;
    n = size(A,1);
    for i = [1:n]
        for j = [1:i-1]
            if A(i,j) ~= 0
                y = 0;
            end
        end
        if A(i,i) == 0
            y = 0;
        end
    end
endfunction

// Función auxiliar para comprobar que la matriz ingresada es triangular inferior
function y = tinf_check(A)
    y = 1;
    n = size(A,1);
    for i = [1:n]
        for j = [i+1:n]
            if A(i,j) ~= 0
                y = 0;
            end
        end
        if A(i,i) == 0
            y = 0;
        end
    end
endfunction

// Resolver un sistema triangular superior Ax = b
//     [a11 a12 ... a1n       [x1       [b1
// A =   0  a22 ... a2n   x =  x2   b =  b2
//      ... ... ... ...        ..        ..
//       0   0  ... ann]       xn]       bn]

// A -> matriz de coeficientes A:nxn
// b -> vector de resultados b:nx1
function x = triang_sup(A,b)
    sz = sz_check(A,b);
    x = %nan;
    if sz == 0 then
        printf("Dimensiones no válidas\n");
    elseif tsup_check(A) == 0
        printf("La matriz debe ser triangular superior\n");
    else
        x(1:sz) = 0;
        for i = flipdim([1:sz],2)
            x(i) = b(i);
            for j = [i+1:sz]
                x(i) = x(i) - A(i,j)*x(j);
            end
            x(i) = x(i)/A(i,i);
        end
    end
endfunction

As = [2,-1,3,-3;0,-1,5,2;0,0,4,1;0,0,0,-3];
bs = [2;-9;-3;9];
xs = triang_sup(As,bs);

// Resolver un sistema triangular inferior Ax = b
//     [a11  0  ...  0        [x1       [b1
// A =  a21 a22 ...  0    x =  x2   b =  b2
//      ... ... ... ...        ..        ..
//      an1 an2 ... ann]       xn]       bn]

// A -> matriz de coeficientes A:nxn
// b -> vector de resultados b:nx1
function x = triang_inf(A,b)
    sz = sz_check(A,b);
    x = %nan;
    if sz == 0 then
        printf("Dimensiones no válidas\n");
    elseif tinf_check(A) == 0
        printf("La matriz debe ser triangular inferior");
    else
        x(1:sz) = 0;
        for i = 1:sz
            x(i) = b(i);
            for j = [1:i-1]
                x(i) = x(i) - A(i,j)*x(j);
            end
            x(i) = x(i)/A(i,i);
        end
    end
endfunction

Ai = [2,0,0;-3,1,0;1,4,-2];
bi = [-2;7;11];
xi = triang_inf(Ai,bi);
