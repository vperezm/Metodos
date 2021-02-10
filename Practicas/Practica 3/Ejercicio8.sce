
// METODO DE PUNTO FIJO

function z = dif(x,y)
    if x > y then
        z = x-y;
    else
        z = y-x;
    end
endfunction

// g/ f(x) = 0 => g(x) = x
function y = m_punto_fijo(g,x0,eps,max_iter)
    iter = 1;
    x = x0;
    while (dif(x,g(x)) > eps) && (iter < max_iter)
        x = g(x);
        iter = iter + 1;
    end
    if iter == max_iter then
        y = %nan;
    else
        y = x;
    end
endfunction

// FUNCIONES

deff("y = g1(x)","y = %e^x/3");
deff("y = g2(x)","y = (%e^x-x)/2");
deff("y = g3(x)","y = log(3*x)");
deff("y = g4(x)","y = %e^x-2*x");

// PRUEBAS
eps = 0.001;
max_iter = 1000;
// g1: x0 in [0,1]
printf("Función: g1(x) = e^x/3\n");
printf("Punto fijo: x = %f\n", m_punto_fijo(g1,0.2,eps,max_iter));
// g2: x0 in [0,1]
printf("Función: g2(x) = (e^x-x)/2\n");
printf("Punto fijo: x = %f\n", m_punto_fijo(g2,0.2,eps,max_iter));
// g3: x0 in [1.1,3]
printf("Función: g3(x) = log(3x)\n");
printf("Punto fijo: x = %f\n", m_punto_fijo(g3,1.2,eps,max_iter));
// g4: x0 in [0.1,1]
printf("Función: g4(x) = e^x-2*x\n");
printf("Punto fijo: x = %f\n", m_punto_fijo(g4,0.2,eps,max_iter));
