clear
// retorna el valor de evaluar la derivada de f de orden n en el punto a
function v = derivar_ci(f,a,h,n)
    if n == 0 then
        v = f(a);
    else
        deff("y = D0F(x)","y = f(x)");
        for i = 1:n
            deff("y = D" +string(i)+ "F(x)", ...
                 "y = (D" +string(i-1)+ "F(x+h) - D" +string(i-1)+ "F(x)) / h");
        end
        deff("y = DnF(x)","y = D" +string(n)+ "F(x)");
        v = DnF(a);
    end
endfunction

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

function y = ex(x)
    y = %e^x;
endfunction

function y = cuad(x)
    y = x^2;
endfunction

a = 2;
h = 0.01;

printf("\nCociente incremental:\n");
printf("f(x)   = \t    e^x \t    x^2\n");
printf("f(2)   = \t %f\t %f\n", derivar_ci(ex, a, h, 0), derivar_ci(cuad, a, h, 0));
printf("f''(2)  = \t %f\t %f\n", derivar_ci(ex, a, h, 1), derivar_ci(cuad, a, h, 1));
printf("f''''(2) = \t %f\t %f\n", derivar_ci(ex, a, h, 2), derivar_ci(cuad, a, h, 2));

printf("\nNumderivative:\n");
printf("f(x)   = \t    e^x \t    x^2\n");
printf("f(2)   = \t %f\t %f\n", derivar_nd(ex, a, h, 0), derivar_nd(cuad, a, h, 0));
printf("f''(2)  = \t %f\t %f\n", derivar_nd(ex, a, h, 1), derivar_nd(cuad, a, h, 1));
printf("f''''(2) = \t %f\t %f\n", derivar_nd(ex, a, h, 2), derivar_nd(cuad, a, h, 2));

// a)
// El error de la implementación con cociente incremental es mayor que el de numderivative.
// b)
// Al ser h un número muy chico, al hacer la resta f(x) - f(x+h) estamos restando dos números muy similares, por lo que puede haber una supresión de cifras significativas.
