function fitness = get_value(position, func_num)
    handle = str2func('cec15_func');
    fitness = feval(handle, position', func_num);
    f = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500];
    fitness = fitness - f(func_num);
end


function fitness = APSO(position, func_num)
    % six benchmark functions used in this paper
    % Unimodal      f1(Sphere)          f2(Quadric)         f3(Rosenbrock)
    % Multimodal    f4(Rastrigrin)      f5(Griewank)        f6(Ackley);
    [np, ~] = size(position);
    fitness = zeros(1, np);
    switch func_num
        case 1
            for index = 1:np
                fitness(index) = Sphere(position(index, :));
            end

        case 2
            for index = 1:np
                fitness(index) = schwefel_p222(position(index, :));
            end

        case 3
            for index = 1:np
                fitness(index) = Quadric(position(index, :));
            end

        case 4
            for index = 1:np
                fitness(index) = Rosenbrockf(position(index, :));
            end

        case 5
            for index = 1:np
                fitness(index) = step(position(index, :));
            end

        case 6
            for index = 1:np
                fitness(index) = Quadric_Noice(position(index, :));
            end
            
        case 7
            for index = 1:np
                fitness(index) = schwefel(position(index, :));
            end

        case 8
            for index = 1:np
                fitness(index) = Rastriginf(position(index, :));
            end

        case 9
            for index = 1:np
                x = position(index, :);
                handle = abs(x) >= 0.5;
                x(handle) = round(2 * x(handle)) / 2;
                fitness(index) = Rastriginf(x);
            end

        case 10
            for index = 1:np
                fitness(index) = Ackley(position(index, :));
            end

        case 11
            for index = 1:np
                fitness(index) = Griewank(position(index, :));
            end
            
        case 12
            for index = 1:np
                fitness(index) = Generalized_Penalized(position(index, :));
            end    
        otherwise 
            disp('func_num if error')
    end
end


%% func_num 1
function value = Sphere(pos)
    value = sum(pos .^ 2);
end

%% func_num 2
function value = Quadric(pos)
    sum1 = 0;
    for index_1 = 1:length(pos)
        sum2 = sum(pos(1:index_1));
        sum1 = sum1 + sum2.^2;
    end
    value = sum1;
end

%% func_num 3
function scores = Rosenbrockf(x)
    scores = 0;
    n = size(x, 2);
    assert(n >= 1, 'Given input X cannot be empty');
    a = 1;
    b = 100;
    for i = 1 : (n-1)
        scores = scores + (b * ((x(:, i+1) - (x(:, i).^2)) .^ 2)) + ((a - x(:, i)) .^ 2);
    end
end

%% func_num 4
function f = Rastriginf(x)
    n = size(x, 2);
    A = 10;
    f = (A * n) + (sum(x .^2 - A * cos(2 * pi * x), 2));
end

%% func_num 5
function value = Griewank(pos)
    sum1 = sum(pos .^ 2);
    prob1 = 1;
    for index = 1:length(pos)
        prob1 = prob1 * cos(pos(index) / sqrt(index));
    end
    value = sum1 / 4000 - prob1 + 1;
end 

%% func_num 6
function scores = Ackley(x)
    n = size(x, 2);
    ninverse = 1 / n;
    sum1 = sum(x .^ 2, 2);
    sum2 = sum(cos(2 * pi * x), 2);
    
    scores = 20 + exp(1) - (20 * exp(-0.2 * sqrt( ninverse * sum1))) - exp( ninverse * sum2);
end

%% func_num 7
% x = [-600,600]
function value = Generalized_Penalized(pos)
    x = pos;
    y = 1 + (x + 1) / 4;
    a = 10; k = 100; m = 4;
    u = zeros(1, length(x));
    handle = x > a;
    u(handle) = k * (x(handle) - a).^m;
    handle = x < -a;
    u(handle) = k * (-x(handle) - a).^m;
    temp = (y(1:end-1) - 1).^2 .* (1 + 10 * sin(pi*y(2:end)).^2);
    sum1 = sum(temp);
    sum2 = sum(u);
    value = (10*sin(pi*y(1))^2 + sum1 + (y(end)-1)^2) * pi / length(x) + sum2;
end

%% func_num 8
% x=[-500, 500] x* = [420.9687, ... , 420.9687]
function scores = schwefel(x)
    n = size(x, 2);
    scores = 418.9829 * n - (sum(x .* sin(sqrt(abs(x))), 2));
end


% x = [-10, 10]
function value = schwefel_p222(x)
    x = abs(x);
    value = sum(x) + prod(x);
end


% x = [-100, 100]
function value = step(x)
    value = sum(floor(x + 0.5) .^ 2);
end

% x = [-1.28, 1.28]
function value = Quadric_Noice(x)
    i = sort(randperm(length(x)));
    value = sum(i .* (x .^ 4)) + rand;
end
    