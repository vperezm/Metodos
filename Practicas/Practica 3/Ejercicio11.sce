// CONSULTAR

clear
// METODO DE NEWTON MULTIVARIABLE
/*
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

function n = f(X)
    x = X(1);
    y = X(2);
    n = 2*x + 3*y^2 + %e^(2*x^2+´y^2);
endfunction

init = [1;1];
eps = 10^(-12);
max_iter = 1000;
printf("\nFunción: f(x,y) = 2*x + 3*y^2 + e^(2*x^2+y^2)\n\n");
res = m_newton_mult(f11,init,eps,max_iter);
printf("Resultado:\tx = %f\ty = %f\n",res(1),res(2));
*/
