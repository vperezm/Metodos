function y = hornerr(p,x0)
    b = flipdim(coeff(p),2);
    for i = [2:length(b)]
        b(i) = b(i) + b(i-1)*x0;
    end
    y = b(length(b));
endfunction

// p(x) = 2 + 3x + 5x^2 + 7x^3
p = poly([2 3 5 7],"x","coeff");
x0 = 2;
// p(x0) = 2 + 3*2 + 5*2^2 + 7*2^3
//       = 2 +  6  +   20  +   56  = 84
b0 = hornerr(p,x0);
b0_sci = horner(p,x0);
printf("p(x) = %s\n", pol2str(p));
printf("Nuestro: p(%d) = %d\n", x0, b0);
printf("Scilab: p(%d) = %d\n", x0, b0_sci);
