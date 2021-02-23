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

printf("Función: (x+1)^(-1)\tIntervalo: [0,1.5]\tSubdivisiones: 10\n");
deff("y = f(x)","y = (x+1)^(-1)");
printf("Método compuesto del trapecio:\n");
yt = trapecioComp(f,0,1.5,10);
disp(yt);
printf("Método compuesto de Simpson:\n");
ys = simpsonComp(f,0,1.5,10);
disp(ys);
printf("Integral:\n");
i = intg(0,1.5,f);
disp(i);
