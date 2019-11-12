    P = 20; D0 = 50; sig = -104; beta = -15; lamF = 320e4; 
    lamD = 5e4; lamA = 1/1800; T = 10000; 
    nbar = pi*power(D0,2)*lamD*(1-exp(-2*0.001* lamA));
    muI_func = @(d2) 2*P*nbar*(power(d2,2)*(-2.5*power(d2,-3.5)) - 2*(d2*8.75*power(d2,-4.5) + 39.375*power(d2,-5.5)))/power(D0,2);
    muI = muI_func(D0) - muI_func(1); %d2 min = 1;

    fun = @(d2) 2*power(P/D0,2)*power(d2,-3);
    temp = integral(fun, 1, D0);            %d2min = 1;
    delI = temp - power(muI/nbar,2);
    delI = delI*sqrt(nbar);
    xi = @(d1) ((P*power(d1,-2.5) - sig*beta)/beta - muI)/delI;
    phi = normcdf(xi,0,1);
    fun2 = @(d1) 2*pi*lamF*d1*phi(d1)*exp(-lamF*pi*power(d1,2));
    out = integral(fun2,1,D0);
