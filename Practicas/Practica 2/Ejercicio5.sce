
function v = derivar_nd(f,a,h,n)
    if n == 0 then
        v = f(a);
    else
        deff("y = D0F(x)","y = f(x)");
        for i = 1:n
            deff("y = D" +string(i)+ "F(x)", ...
                 "y = numderivative(D" +string(i-1)+ "F,x,h)");
        end
        deff("y = DnF(x)","y = D" +string(n)+ "F(x)");
        v = DnF(a);
    end
endfunction

// Calcula el valor del polinomio de Taylor de grado n de una función f en un punto v
function y = taylor(f,n,a,v)
    coeficientes(1) = 0;
    for i = 1:n
        coeficientes(i) = (1/factorial(i-1))*derivar_nd(f,a,0.001,i-1);
    end
    y = horner(poly([coeficientes],"x","coeff"),(v-a));
endfunction

function y = ex(x)
    y = %e^x;
endfunction

function y = cuad(x)
    y = x^2;
endfunction

printf("Función: e^x\tValor: x = 1\n");
printf("Esperado: e^1 = %f\n", %e^1);
printf("Talylor: e^1 = %f\n", taylor(ex,7,0,1));
