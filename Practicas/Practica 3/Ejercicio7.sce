clear
// a)
function z = dif(x,y)
    if x > y then
        z = x-y;
    else
        z = y-x;
    end
endfunction

// METODO DE PUNTO FIJO
// f/ g(x) = 0 => f(x) = x
function y = m_punto_fijo(f,x0,eps,max_iter)
    iter = 1;
    x = x0;
    while (dif(x,f(x)) > eps) && (iter < max_iter)
        x = f(x);
        iter = iter + 1;
    end
    if iter == max_iter then
        y = %nan;
    else
        y = x;
    end
endfunction

// w^2 = g*d*tanh(hd)
// => d = w^2/(g*tanh(hd))
// w = 2*%pi/T
// d = 2*%pi/l
// g = 9.8 m/s^2
// h = 4 m
// T = 5 s
function d = f1(x)
    T = 5;
    w = 2*%pi/T;
    g = 9.8;
    h = 4;
    d = w^2/(g*tanh(h*x));
endfunction

// El método de punto fijo dará la solución en relación de d
// Debo despejar l
// d = 2*%pi/l
// => l = 2*%pi/d
function l = d2l(d)
    l = 2*%pi/d;
endfunction

max_iter = 1000;
res_a = m_punto_fijo(f1,1,0.1,max_iter);
printf("Longitud de onda\n");
printf("Punto fijo: d = %f\tl = %f\n", res_a, d2l(res_a));

// b)

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

// w^2 = g*d*tanh(hd)
// => 0 = w^2 - g*d*tanh(hd)
function y = f2(x)
    T = 5;
    w = 2*%pi/T;
    g = 9.8;
    h = 4;
    y = w^2-(g*x*tanh(h*x));
endfunction

res_b = m_newton(f2,res_a,0.0001,max_iter);
printf("Newton: d = %f\tl = %f\n", res_b, d2l(res_b));
