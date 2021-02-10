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

// ----------------------------------------------------------------------------//
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

// ----------------------------------------------------------------------------//
// METODO DE NEWTON
// f y f' continuas en [a,b]
// alfa / f(alfa) = 0
// x0 en [a,b] cercano a alfa y f'(x0) != 0
function y = m_newton(f,x0,eps,max_iter)
    iter = 1;
    x = x0;
    while (abs(f(x)) > eps) && (iter < max_iter)
        x = x - f(x) / numderivative(f,x);
        iter = iter + 1;
    end
    if iter == max_iter then
        y = %nan;
    else
        y = x;
    end
endfunction

// ----------------------------------------------------------------------------//
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

// ----------------------------------------------------------------------------//
// METODO DE NEWTON MULTIVARIABLE

function y = m_newton_mult(fn,X,eps,max_iter)
    iter = 1;
    ant = X;
    act = X;
    while (norm(act-ant) > eps | iter == 1) && (iter < max_iter)
        ant = act;
        J = numderivative(fn,act);
        act = act - inv(J)*fn(act);
        iter = iter + 1;
    end
    if iter == max_iter then
        y = %nan;
    else
        y = act;
    end
endfunction
