function [bestfit] = CMA_ES(FUNC, DIM, UB, LB)
%% parameter settings
sigma = 1.0;
chin = DIM^0.5*(1 - 1/(4*DIM) + 1/(21*DIM^2));
cs = 4/(DIM + 4);
cc = 4/(DIM + 4);
ccov = 2/(DIM+2^0.5)^2;
damp = 1/cs + 1;

lambda = 4+floor(3*log(DIM));
mu = floor(lambda/2);
weights = log((lambda+1)/2) - log(1:mu)';

NP = lambda;

B =  eye(DIM); D = eye(DIM); C = B*D*(B*D)';
pc = zeros(DIM,1); ps = zeros(DIM,1);
cw = sum(weights)/norm(weights);

x = zeros(DIM, NP);
z = zeros(DIM, NP);
fitness = zeros(NP, 1);
fes = 0;


for k=1:lambda
    % repeat the next two lines until arx(:,k) is feasible
    x(:,k) = rand(1, DIM)*(UB - LB) + LB;
    fitness(k) = benchmark_2015(x(:,k)', FUNC);
    fes = fes+1;
end

[fitness, index] = sort(fitness); % minimization
xmeanw = x(:,index(1:mu))*weights/sum(weights);

endfit = 10^-8;
maxfes = DIM*10000;


bestfit = 10^10;



%% main loop

while fes < maxfes && bestfit > endfit
    for k=1:lambda
        % repeat the next two lines until arx(:,k) is feasible
        z(:,k) = randn(DIM,1);
        x(:,k) = xmeanw + sigma * (B * D * z(:,k)); % Eq.(13)
        x(:, k) = min(max(x(:, k), LB), UB);
        fitness(k) = benchmark_2015(x(:,k)', FUNC);
        fes = fes+1;
    end
    
    [fitness, index] = sort(fitness); % minimization
    xmeanw = x(:,index(1:mu))*weights/sum(weights);
    zmeanw = z(:,index(1:mu))*weights/sum(weights);
    
    pc = (1 - cc)*pc + (sqrt(cc*(2-cc))*cw)*(B*D*zmeanw);
    C = (1 - ccov)*C + ccov*pc*transpose(pc);
    
    ps = (1 - cs)*ps + (sqrt(cs*(2-cs))*cw)*(B*zmeanw);
    sigma = sigma*exp((norm(ps) - chin)/(chin*damp));
    
    [B, D] = eig(C);
    D = diag(sqrt(diag(D)));
    bestfit = min(fitness);
    
    disp(['step_search:', num2str(fes),  ' best fit: ', num2str(bestfit)]);
    
end

end


