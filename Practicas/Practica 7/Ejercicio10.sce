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


// Chebyshev(n) calcula las raíces de un polinomio de Chebyshev para n nodos interpolantes
function rot = Chebyshev(n)
    for i = 1:n
        rot(i) = cos(((2*i)-1)*%pi/(2*n));
    end
endfunction

function [x,P] = pol_int_cheb(n)
    x = Chebyshev(n+1);
    x = x';
    for i = 1:n+1
        y(1,i) = ej10(x(i));
    end
    P = polydifdiv(x,y);
    
endfunction

function y = ej10(x)
    y = %e.^x;
endfunction

// a)
[x,P] = pol_int_cheb(3);
disp(P);

// b)
function y = err_ej10(x,P)
    n = length(x);
    for i = 1:n
        y(1,i) = abs(ej10(x(i))-horner(P,x(i)));
    end
endfunction

x = [-1:0.01:1];
y = err_ej10(x,P);
plot(x,y);
