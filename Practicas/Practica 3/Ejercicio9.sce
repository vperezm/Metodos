// METODO DE NEWTON MULTIVARIABLE
// (cantidad fija de iteraciones)

// fn -> sistema de n ecuaciones no lineales
// Xn -> n variables
// N  -> numero de iteraciones

function y = m_newton_mult(fn,X,N)
    Xn = X;
    for i = 1:N
        J = numderivative(fn,Xn);
        y = Xn - inv(J)*fn(Xn);
        Xn = y;
    end
endfunction

function n = f9(X)
    x = X(1);
    y = X(2);
    n = [1 + x^2 - y^2 + %e^x*cos(y); 2*x*y + %e^x*sin(y)];
endfunction

init = [-1;4];
iter = 5;
res = m_newton_mult(f9,init,iter);
printf("Sistema de ecuaciones:\n1 + x^2 - y^2 + e^x*cos(y)\n2*x*y + e^x*sin(y)\n\n");
printf("Resultados:\nx = %f\ny = %f\n",res(1),res(2));
