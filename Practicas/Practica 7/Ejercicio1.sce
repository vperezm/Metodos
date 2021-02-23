clear

// METODO DE INTERPOLACION DE LAGRANGE
// P -> conjunto de pares de elementos (x,f(x))
// x -> punto a aproximar

function y = interpolacion_lagrange(P,x)
    n = size(P,"r");
    y = 0;
    
    for k = 1:n
        // Calculamos la función Lk(x)
        Lk = 1;
        for i = 1:n
            if i <> k
                Lk = Lk * (x - P(i,1)) / (P(k,1) - P(i,1));
            end
        end
        // Sumamos el término k del polinomio
        y = y + Lk*P(k,2);
    end
endfunction

// METODO DE LAS DIFERENCIAS DIVIDIDAS DE NEWTON

// Conjunto de diferencias divididas de Newton
// P -> conjunto de pares de elemntos (x,f(x))
function D = DDs(P)
    n = size(P,"r")
    for i = 1:n
        D(i,1) = P(i,2)
    end
    
    for i = 2:n
        for j = 1:n-(i-1)
            D(j,i) = ( D(j+1,i-1) - D(j,i-1)) / ( P(j+i-1,1) - P(j,1) )
        end
    end
    //disp(D)
endfunction

// Interpolación de Newton
function y = interpolacion_newton(P,x)
    n = size(P,"r");
    y = P(1,2);
    D = DDs(P);
    for i = 2:n
        pr = 1;
        for j = 1:i-1
            pr = (x - P(j,1)) * pr;
        end
        y = y + pr*D(1,i);
    end
endfunction

//------------------------------------------------------------------------------//

x = 1/3;
y = 1.395612425;
P1 = [0,1.0;0.6,1.8221];
P2 = [0,1.0;0.2,1.2214;0.4,1.4918;0.6,1.8221];

printf("INTERPOLACIÓN DE LAGRANGE:\n");
yl1 = interpolacion_lagrange(P1,x);
printf("Lineal: %f\tError: %f\n",yl1,abs(y-yl1));
yl2 = interpolacion_lagrange(P2,x);
printf("Cúbica: %f\tError: %f\n",yl2,abs(y-yl2));

printf("-----------------------------\n");

printf("INTERPOLACIÓN DE NEWTON:\n");
yn1 = interpolacion_newton(P1,x);
printf("Lineal: %f\tError: %f\n",yn1,abs(y-yn1));
yn2 = interpolacion_newton(P2,x);
printf("Cúbica: %f\tError: %f\n",yn2,abs(y-yn2));

printf("-----------------------------\n");

// Cálculo de las cotas del error
// Sea f(x) = e^x. Luego f'(x) = e^x.
// También, f es estrictamente creciente en todo su dominio.
// Considerando el intervalo [0,0.6], min(f(x)) = f(0) y max(f(x)) = f(0.6)
function [mi,ma] = error_ej1(P,x)
    n = size(P,1);
    mi = P(1,2);
    ma = P(n,2);
    mul = 1;
    for i = 1:n
        mul = mul * (x - P(i,1));
    end
    mul = mul / factorial(n+1);
    mi = abs(mi*mul);
    ma = abs(ma*mul);
endfunction

printf("Cotas del error:\n");
[mi1,ma1] = error_ej1(P1,x);
printf("Lineal: %f <= E <= %f\n",mi1,ma1);
[mi2,ma2] = error_ej1(P2,x);
printf("Cúbica: %f <= E <= %f\n",mi2,ma2);
