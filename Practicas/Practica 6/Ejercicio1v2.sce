
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

//------------------------------------------------------------------------------//

// a)

A = [1,0,0;-1,0,1;-1,-1,2];

printf("--------------------------------------------\n");
printf("Matriz:\n");
disp(A);
printf("Autovalores:\n");
disp(spec(A));
gershgorin(A);

// b)

B = [1,0,0;-0.1,0,0.1;-0.1,-0.1,2];

printf("--------------------------------------------\n");
printf("Matriz:\n");
disp(B);
printf("Autovalores:\n");
disp(spec(B));
gershgorin(B);

// c)

C = [1,0,0;-0.25,0,0.25;-0.25,-0.25,2];

printf("--------------------------------------------\n");
printf("Matriz:\n");
disp(C);
printf("\nAutovalores:\n");
disp(spec(C));
gershgorin(C);


// d)

D = [4,-1,0;-1,4,-1;-1,-1,4];

printf("--------------------------------------------\n");
printf("Matriz:\n");
disp(D);
printf("Autovalores:\n");
disp(spec(D));
gershgorin(D);

// e)

E = [3,2,1;2,3,0;1,0,3];

printf("--------------------------------------------\n");
printf("Matriz:\n");
disp(E);
printf("Autovalores:\n");
disp(spec(E));
gershgorin(E);


// f)

F = [4.75,2.25,-0.25;2.25,4.75,1.25;0.25,1.25,4.75];

printf("--------------------------------------------\n");
printf("Matriz:\n");
disp(F);
printf("Autovalores:\n");
disp(spec(F));
gershgorin(F);
