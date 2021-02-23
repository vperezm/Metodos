clear
clc

// FactoRQ(A) Factoriza una matriz A en Q,R (usado en mínimos cuadrados)
function [Q,R] = FactoRQ(A)   // A debe tener columnas LI
    sz = size(A,2)
    Q(:,1) = A(:,1)/norm(A(:,1))
    V(1) = norm(A(:,1))
    for i = 2:sz
        suma = 0
        for j = 1:i-1
            suma = suma + (A(:,i)'*Q(:,j))*Q(:,j)
        end
        Q(:,i) = A(:,i) - suma
        V(i) = norm(Q(:,i))
        Q(:,i) = Q(:,i)/V(i)  
    end
    
    R = diag(V)
    
    for i = 1:sz
        for j = i+1:sz
           R(i,j) = A(:,j)'*Q(:,i) 
        end
    end
endfunction


//minimoscuad(xi,y,gr): Toma un conjunto de puntos xi, y un conjunto de f(xi)=y, y un grado de polinomio. Aproxima un polinomio por minimos cuadrados y devuelve el error.
function [p, err] = minimoscuad(xi, y, gr)
    szy = size(y,1);
    n = gr+1;
    A = eye(szy, n);
    
    for j = 1:n
        for i = 1:szy
            A(i, j) = xi(i)**(j-1);    
        end
    end
    
    [Q,R] = FactoRQ(A);
    
    b = Q'*y
    
    sz = size(R,1)
    
    x(sz) = b(sz)/R(sz,sz)
    
    for i = 1:sz-1 //ResUELVE QR
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*R(sz-i,sz-j+1)   
        end
        x(sz-i) = (b(sz-i)-suma)/R(sz-i,sz-i) 
    end
    
    p = poly(x, 'x', "coeff")
    
    E = A*x-y;
    err = E'*E;
    
endfunction

x = [0;0.15;0.31;0.5;0.6;0.75];
y = [1;1.004;1.31;1.117;1.223;1.422];

printf("----------------------------------\n");
[p1,e1] = minimoscuad(x,y,1);
printf("Grado 1:\n\n");
printf("Polinomio:\n");
disp(p1);
printf("Error:\n");
disp(e1);

printf("----------------------------------\n");
[p2,e2] = minimoscuad(x,y,2);
printf("Grado 2:\n\n");
printf("Polinomio:\n");
disp(p2);
printf("Error:\n");
disp(e2);

printf("----------------------------------\n");
[p3,e3] = minimoscuad(x,y,3);
printf("Grado 3:\n\n");
printf("Polinomio:\n");
disp(p3);
printf("Error:\n");
disp(e3);
