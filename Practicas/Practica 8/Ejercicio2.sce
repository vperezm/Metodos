clear
clc


// Método del trapecio compuesto para n subdivisiones
function y = trapecioComp(f,a,b,n)
    h = (b-a)/n;
    y = f(a)*1/2 + f(b)*1/2;
    for j = 1:n-1
        y = y + f(a+j*h);
    end
    y = y*h;
endfunction

// Función auxiliar para imprimir los resultados
function x = display(f,a,b,n)
    printf("Método compuesto del trapecio:\n");
    yt = trapecioComp(f,a,b,n);
    disp(yt);
    printf("Integral:\n");
    i = intg(a,b,f);
    disp(i);

    x = 0;
endfunction

printf("\nb)\n");
printf("Función: x^3\tIntervalo: [0,2]\tSubdivisiones: 4\n");
deff("y = fb(x)","y = x^3");
x = display(fb,0,2,4);

printf("\nc)\n");
printf("Función: x*(1+x^2)^(1/2)\tIntervalo: [0,3]\tSubdivisiones: 6\n");
deff("y = fc(x)","y = x*(1+x^2)^(1/2)");
x = display(fc,0,3,6);

printf("\nd)\n");
printf("Función: sin(pi*x)\tIntervalo: [0,1]\tSubdivisiones: 8\n");
deff("y = fd(x)","y = sin(%pi*x)");
x = display(fd,0,1,8);

printf("\ne)\n");
printf("Función: x*sin(x)\tIntervalo: [0,2*pi]\tSubdivisiones: 8\n");
deff("y = fe(x)","y = x*sin(x)");
x = display(fe,0,2*%pi,8);

printf("\nf)\n");
printf("Función: x^2*e^x\tIntervalo: [0,1]\tSubdivisiones: 8\n");
deff("y = ff(x)","y = (x^2)*(%e^x)");
x = display(ff,0,1,8);
