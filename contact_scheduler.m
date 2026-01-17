function schedule = contact_scheduler(dict, dt)
    run_time = 10; % 10 seconds
    double_contact_time = .25;
    mode = dict('movement_type').value;
    switch mode
        case "standing"
            schedule = ones(run_time/dt, 2); % Standing schedule
        case "walking"
            right = repmat([0,1], 1/dt , 1); % Right foot on ground for 1s
            left = repmat([1,0], 1/dt , 1); % Left foot on ground for 1s
            standing = repmat([1,1], double_contact_time/dt, 1); % Both feet contacted for .25s
            schedule = repmat([standing; right; standing; left;], run_time / 2.5, 1);
        case "running"
    end
end
