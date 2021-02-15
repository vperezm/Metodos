// Cálculo del determinante mediante el método de Gauss

function y = det_gauss(A)
    [nA,mA] = size(A);
    n = nA;
    a = A;
    for k=1:n
        for i=k+1:n
            a(i,k+1:n) = a(i,k+1:n) - a(k,k+1:n)*a(i,k)/a(k,k);
            a(i,1:k) = 0;  // no hace falta para calcular la solución x
        end;
    end;
    disp(a);
    y = 1;
    for i = 1:n
        y = y * a(i,i);
    end
endfunction

A = [1,2,4,1;2,5,7,-1;-1,1,-1,2;3,5,9,3];
y = det_gauss(A);
printf("det(A) = %f\n",y);
printf("Scilab: %f\n", det(A));
