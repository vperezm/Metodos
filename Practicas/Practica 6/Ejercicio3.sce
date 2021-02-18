//gershgorin(A). Muestra las cotas de los autovalores de A dada dicha matriz.

function gershgorin(A)
    sz = size(A, 1);
    
    for i = 1:sz
        suma = 0;
        for j = 1:sz
            if (i <> j)
                suma = suma + abs(A(i,j));
            end;
        end;
            mprintf("|lambda - %f| <= %f\n", A(i,i), suma);
    end;
endfunction


// poly([A], "x") -> polinomio caracteristico de la matriz A 
//(det(lambda*I - A) = p(lambda))

function ej_3(A)
   sz = size(A, 1);
   for k = 0:10
       mprintf("k = %d\n", k)
       
       A(sz,sz) = 1 + 0.1*k;
       p = poly([A], "x");
       x = roots(p);
       disp(x)
       av = spec(A);
       disp(av)
       gershgorin(A);
   end
endfunction

A = [1,-1,0;-2,4,-2;0,-1,1];
// el ultimo elemento corresponde a A(3,3), que será modificado en la función para representar los distintos valores de 1 + eps
ej_3(A);
