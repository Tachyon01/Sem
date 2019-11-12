clc;
clear all;
close all;

%Inputs
r = 150;%input('Radius of considered area');
D0 = 50;%input('Radius of interference area D0');
P = 20;%input('TN transmit power P');
% gd = input('Path loss model g(d)');
T = 10000;%input('Total time T');
L = 256;%input('Transaction packet lenght L');
lamA = 1/1800;%input('Transaction arrival density lambdaA');
sig = -104;%input('Noise power sigma');


%initialization
lamF=320e4;
beta = -15;

%initializations done
i=5; out = zeros(1,14);
for lamD = 5e4:1e4:14e4
    
    nbar = pi*power(D0,2)*lamD*(1-exp(-2*0.001* lamA));
    muI_func = @(d2) 2*P*nbar*(power(d2,2)*(-2.5*power(d2,-3.5)) - 2*(d2*8.75*power(d2,-4.5) + 39.375*power(d2,-5.5)))/power(D0,2);
    muI = muI_func(D0) - muI_func(1); %d2 min = 1;

    fun = @(d2) 2*power(P/D0,2)*power(d2,-3);
    temp = integral(fun, 1, D0);            %d2min = 1;
    delI = temp - power(muI/nbar,2);
    delI = delI*sqrt(nbar);

    xi = @(d1) ((P*power(d1,-2.5) - sig*beta)/beta - muI)/delI;					%beta is different for different runs
%     N0 = pi*r*r*lamD*10e4;
%     muM = N0*lamA*T;
    phi = normcdf(xi,0,1);
    fun2 = @(d1) 2*pi*lamF*d1*phi(d1)*exp(-lamF*pi*power(d1,2));
    out = integral(fun2,1,D0);
	i=i+1;
end
plot(out);
