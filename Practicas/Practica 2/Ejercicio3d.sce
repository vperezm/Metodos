clear

function y = horner_der(p,x0) // y = [p(x0) p'(x0)]
    b = flipdim(coeff(p),2);
    n = length(b)-1;
    q = 0;
    for i = [1:length(b)]
        if i > 1
            b(i) = b(i) + b(i-1)*x0;
        end
        if i < length(b)
            q = q + b(i)*x0^(n-i) 
        end
    end
    y = [b(length(b)) q];
endfunction

// p(x) = 2 + 3x + 5x^2 + 7x^3
p = poly([2 3 5 7],"x","coeff");
x0 = 2;
// p(x0) = 2 + 3*2 + 5*2^2 + 7*2^3
//       = 2 +  6  +   20  +   56  = 84
// p'(x) = 3 + 10x + 21x^2
// p'(x0) = 3 + 20 + 84 = 107
t = horner_der(p,x0);
printf("p(x) = %s\n", pol2str(p));
printf("p(%d) = %d  (expected 84)\n", x0, t(1));
printf("p''(%d) = %d  (expected 107)\n", x0, t(2));
