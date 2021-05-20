function f = benchmark_2015(x, func_num)
    fhd = str2func('cec15_func');
    f = feval(fhd, x', func_num) - bias(func_num);
    f = f';
end

function outcome = bias(func_num)
    outcome = func_num * 100;
end
  