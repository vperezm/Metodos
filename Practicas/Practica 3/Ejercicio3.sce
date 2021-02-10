// METODO DE LA SECANTE

function y = m_secante(f,a,b,eps,max_iter)
    iter = 1;
    actual = b - f(b)*(b-a)/(f(b)-f(a));
    anterior = b;
    while (abs(f(actual)) >= eps) && (iter < max_iter)
        siguiente = actual - f(actual)*(actual-anterior)/(f(actual)-f(anterior));
        anterior = actual;
        actual = siguiente;
        iter = iter + 1;
    end
    if iter == max_iter then
        y = %nan;
    else
        y = actual;
    end
endfunction

function y = f(x)
    y = x^2/4 - sin(x);
endfunction

eps = 0.01;
max_iter = 1000;
printf("Función: f(x) = x^2/4 - sin(x)\n");
printf("Raíz 1: x = %f\n", m_secante(f,1,4,eps,max_iter));
