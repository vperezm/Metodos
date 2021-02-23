clear
clc

// Regla del trapecio simple extendida
// Funciones de dos variables f(x,y)
// Intervalos: x in [a,b], y in [c,d]
function y = trapecioExt(f,a,b,c,d)
    hy = d-c;
    hx = b-a;
    y = (f(c,a)+f(c,b)+f(d,a)+f(d,b))*hx*hy/4;
endfunction

printf("Función: sin(x+y)\tIntervalos: x in [0,1], y in [0,2]\n");
deff("y = f(x,y)","y = sin(x+y)");
y = trapecioExt(f,0,1,0,2);
disp(y);
