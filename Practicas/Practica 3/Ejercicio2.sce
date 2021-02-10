// METODO DE LA BISECCION
// f continua en [a,b] tal que f(a)f(b) < 0
// eps tolerancia del error
// max_iter maximo de iteraciones

function y= m_biseccion(f,a,b,eps,max_iter)
    iter = 1;
    c = (a+b)/2;
    while (b-c > eps) && (f(c) ~= 0) && (iter < max_iter) 
        if f(a)*f(c) < 0
            b = c;
        else
            a = c;
        end
        c = (a+b)/2;
        iter = iter + 1;
    end
    if iter == max_iter then
        y = %nan;
    else
        y = c;
    end
endfunction

function y = f1(x)
    y = sin(x)-x^2/2;
endfunction

function y = f2(x)
    y = x^4 - %e^(-x);
endfunction

function y = f3(x)
    y = x - log(x) - 1;
endfunction

eps = 0.1;
max_iter = 1000;
printf("Función: f1(x) = sin(x)-x^2/2\n");
printf("Raíz 1: x = %f\n", m_biseccion(f1,-1,1,eps,max_iter));
printf("Raíz 2: x = %f\n", m_biseccion(f1,1,3,eps,max_iter));
printf("Función: f2(x) = x^4 - e^(-x)\n");
printf("Raíz 1: x = %f\n", m_biseccion(f2,-10,-7,eps,max_iter));
printf("Raíz 2: x = %f\n", m_biseccion(f2,-2,0,eps,max_iter));
printf("Raíz 3: x = %f\n", m_biseccion(f2,0,2,eps,max_iter));
printf("Función: f3(x) = x - log(x) - 1\n");
printf("Raíz 1: x = %f\n", m_biseccion(f3,0.1,0.5,eps,max_iter));
printf("Raíz 1: x = %f\n", m_biseccion(f3,0.5,1.5,eps,max_iter));

