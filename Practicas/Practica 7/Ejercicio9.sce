xdel
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

function y = ej9(x)
    y = 1/(1+x.^2);
endfunction

function y = values(x)
    n = length(x);
    for i = 1:n
        y(i) = ej9(x(i));
    end
    y = y';
endfunction

function x = aux(n)
    x = [-5:(10/n):5];
endfunction

function y = err_ej9(x,P)
    n = length(x);
    for i = 1:n
        y(1,i) = abs(ej9(x(i))-horner(P,x(i)));
    end
endfunction

function P = forplot(n)
    x = aux(n);
    y = values(x);
    P = polydifdiv(x,y);
endfunction

/*
// Grado 2
P = forplot(2);
x = [-5:0.01:5];
plot(x,err_ej9(x,P));
*/
/*
// Grado 4
P = forplot(4);
x = [-5:0.01:5];
plot(x,err_ej9(x,P));
*/
/*
// Grado 6
P = forplot(6);
x = [-5:0.01:5];
plot(x,err_ej9(x,P));
*/
/*
// Grado 10
P = forplot(10);
x = [-5:0.01:5];
plot(x,err_ej9(x,P));
*/

// Grado 14
P = forplot(14);
x = [-5:0.01:5];
plot(x,err_ej9(x,P));

