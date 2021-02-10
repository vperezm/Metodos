clear
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

function n = f10(X)
    x = X(1);
    y = X(2);
    n = [x^2 + x*y^3 - 9; 3*x^2*y - 4 - y^3];
endfunction

X1 = [1.2;2.5];
X2 = [-2;2.5];
X3 = [-1.2;-2.5];
X4 = [2;-2.5];
eps = 0.000001;
max_iter = 1000;

printf("\nSistema de ecuaciones:\n");
printf("x^2 + x*y^3 - 9\n3*x^2*y - 4 - y^3\n\n");
printf("a)\n");
res1 = m_newton_mult(f10,X1,eps,max_iter);
check1 = f10(res1);
printf("Valores iniciales:\tx = %f\ty = %f\t",X1(1),X1(2));
printf("Resultado:\tx = %f\ty = %f\n",res1(1),res1(2));
printf("Check: x^2 + x*y^3 - 9 = %f\t3*x^2*y - 4 - y^3 = %f\n\n",check1(1),check1(2));
printf("b)\n");
res2 = m_newton_mult(f10,X2,eps,max_iter);
check2 = f10(res2);
printf("Valores iniciales:\tx = %f\ty = %f\t",X2(1),X2(2));
printf("Resultado:\tx = %f\ty = %f\n",res2(1),res2(2));
printf("Check: x^2 + x*y^3 - 9 = %f\t3*x^2*y - 4 - y^3 = %f\n\n",check2(1),check2(2));
printf("c)\n");
res3 = m_newton_mult(f10,X3,eps,max_iter);
check3 = f10(res3);
printf("Valores iniciales:\tx = %f\ty = %f\t",X3(1),X3(2));
printf("Resultado:\tx = %f\ty = %f\n",res3(1),res3(2));
printf("Check: x^2 + x*y^3 - 9 = %f\t3*x^2*y - 4 - y^3 = %f\n\n",check3(1),check3(2));
printf("d)\n");
res4 = m_newton_mult(f10,X4,eps,max_iter);
check4 = f10(res4);
printf("Valores iniciales:\tx = %f\ty = %f\t",X4(1),X4(2));
printf("Resultado:\tx = %f\ty = %f\n",res4(1),res4(2));
printf("Check: x^2 + x*y^3 - 9 = %f\t3*x^2*y - 4 - y^3 = %f\n\n",check4(1),check4(2));
