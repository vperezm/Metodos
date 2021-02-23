clear
clc

// Método del trapecio simple
function y = trapecio(f,x0,x1)
    y = (x1-x0)*(f(x0)+f(x1))/2;
endfunction

// Método del trapecio compuesto para n subdivisiones
function y = trapecioComp(f,a,b,n)
    h = (b-a)/n;
    y = f(a)*1/2 + f(b)*1/2;
    for j = 1:n-1
        y = y + f(a+j*h);
    end
    y = y*h;
endfunction

// Método de Simpson
function y = simpson(f,a,b)
    h = (b-a)/2;
    y = (f(a) + 4*f(a+h) + f(b))*h/3;
endfunction

// Método compuesto de Simpson para n subdivisiones
function y = simpsonComp(f,a,b,n)
    h = (b-a)/n;
    y = f(a) + f(b);
    aux = 4;
    for j = 1:n-1
        y = y + aux*f(a+j*h);
        if aux == 4
            aux = 2;
        else
            aux = 4;
        end
    end
    y = y*h/3;
endfunction


// Función auxiliar para imprimir los resultados
function x = display(f,a,b)
    printf("Método del trapecio:\n");
    yt = trapecio(f,a,b);
    disp(yt);
    printf("Método de Simpson:\n");
    ys = simpson(f,a,b);
    disp(ys);
    printf("Integral:\n");
    i = intg(a,b,f);
    disp(i);

    x = 0;
endfunction

printf("\na)\n");
printf("Función: ln(x)\tIntervalo: [1,2]\n");
x = display(log,1,2);

printf("\nb)\n");
printf("Función: x^(1/3)\tIntervalo: [0,0.1]\n");
deff("y = fb(x)","y = x^(1/3)");
x = display(fb,0,0.1);

printf("\nc)\n");
printf("Función: sin^2(x)\tIntervalo: [0,pi/3]\n");
deff("y = fc(x)","y = (sin(x))^2");
x = display(fc,0,%pi/3);

// intg(a,b,f) calcula la integral definida
// a -> extremo inferior de la integral
// b -> extremo superior de la integral
// f -> la función a ser integrada
