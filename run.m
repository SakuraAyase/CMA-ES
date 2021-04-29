FUNC_LIST = 15;
RUN_TIMES = 5;
record_list= cell(1, 4);
dim_list = [10, 30, 50, 100];



for i = 1:4
    DIM = dim_list(i);

    record = zeros(FUNC_LIST, RUN_TIMES);
    % parpool('local')
    for func_num = 1:FUNC_LIST
        rng(1);
        for runtime = 1:RUN_TIMES
            
            record(func_num, runtime) = CMA_ES(func_num, DIM, 100, -100);
        end
    end
    record_list{i} = record;
end
