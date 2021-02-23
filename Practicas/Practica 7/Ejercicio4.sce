clear
clc

// Diferencia dividida de orden k
// x,y vectores de valores donde (x(i),y(i)) son los datos para la interpolación y D es la diferencia dividida de orden k correspondiente
function D = difdivk(x,y)
    k = length(x);
    if k == 2 then
        D = (y(2)-y(1))/(x(2)-x(1));
    else
        D = (difdivk(x(2:k),y(2:k))-difdivk(x(1:k-1),y(1:k-1)))/(x(k)-x(1));
    end
endfunction

// Polinomio de interpolación por diferencias divididas
// v,w vectores de valores donde (v(i),w(i)) son los datos dados para la interpolación. P es el polinomio de interpolación por diferencias divididas correspondiente
function P = polydifdiv(v,w)
    n = length(v);
    P = w(1);
    for k = 2:n
        P = P + difdivk(v(1:k),w(1:k))*poly(v(1:k-1),"x",["roots"]);
    end
endfunction

//------------------------------------------------------------------------------//

// Definimos los vectores x e y correspondientes a los valores dados en la tabla:
x = [2,2.1,2.2,2.3,2.4,2.5];
y = [0.2239,0.1666,0.1104,0.0555,0.0025,-0.0484];

// Calculamos el polinomio de interpolacion por diferencias de newton usando los puntos dados:
P = polydifdiv(x,y);

// Evaluamos dicho polinomio en x = 2.15 y en x = 2.35
aprox1 = horner(P,2.15);
aprox2 = horner(P,2.35);

printf("\nJ0(2.15) = %f\n", aprox1);
printf("J0(2.35) = %f\n", aprox2);
