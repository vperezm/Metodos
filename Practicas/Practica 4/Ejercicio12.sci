function [U, ind] = Cholesky(A)
eps = 1.0e-8
n = size(A,1)
U = zeros(n,n)
for k = 1:n
    if k==1 then
            t = A(k,k)
    else 
            t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    end

    if t <= eps
        printf("Matriz no definida positiva.\n")
        ind = 0
        return
    end
    U(k,k)= sqrt(t)
    for j = k+1:n
        if k==1 then 
                    U(k,j) = A(k,j)/U(k,k)
        else 
                    U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
        end
    end
end
ind = 1
endfunction

// Sistemas triangulares superiores
function x = sust_reg(A,b)
    n = size(A,1);
    x(n) = b(n)/A(n,n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i,i+1:n)*x(i+1:n))/A(i,i);
    end
endfunction

// Sistemas triangulares inferiores
function x = sust_prog(A,b)
    n = size(A,1);
    x(1) = b(1)/A(1,1);
    for i = 2:n
        x(i) = (b(i) - A(i,1:i-1)*x(1:i-1))/A(i,i);
    end
endfunction

function x = chol_solver(A,b)
    [R,ind] = Cholesky(A);
    if ind == 0 then
        error('No se puede calcular');
        abort;
    end
    g = sust_prog(R',b);
    x = sust_reg(R,g);
endfunction

A = [16,-12,8;-12,18,-6;8,-6,8];
b = [76;-66;46];
x = chol_solver(A,b);
disp(x);
printf("Verificación:\n");
printf("Ax:\n");
disp(A*x);
printf("b:\n");
disp(b);
