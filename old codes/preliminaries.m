% MCMCMC method and its TMCMC analog

fplot(@(x) pdfmix(x,0,10,1,1,0.5),[-4 12])
xlabel('x');
ylabel('density');

handle=@(x) pdfmix(x,0,8,1,1,0.5);
beta=1;
fplot(@(x) power(handle(x),beta),[-4 12]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1=repmat(0,1,d);
mu2=repmat(4,1,d);
Sigma1=eye(d);
Sigma2=eye(d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting the sampled cold chain after MC3

t=length(Out);
plot(1:length(Out),Out);
ylabel('simulated value');
xlabel('iterate');
