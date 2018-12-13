function identify_INT_vars(vars, vals)
    
    #count int variables
    c = 0
    int_vars = String[]
    int_vals = Float64[]
    
    for (i,var) in enumerate(vars)
        key = split(var,"_")[1]
        if key == "INT"
            push!(int_vars,var)
            push!(int_vals,vals[i])
        end
    end
        return int_vars,int_vals
    end

function nearest_ints(val)
    rounded = round(val)
    tol = 1e-15
    if rounded >= val
        return (rounded, rounded-1+tol)
    else
        return (rounded+1,rounded)
    end
end

function cleanup_input(varname::String)
    trim = split(varname, '.')[1]
    elem = split(trim, '_')
    s = ""
    #@show varname, elem
    for i in 1:length(elem)-1
        s = s*string(elem[i])*'_'
    end
    s = s[1:end-1]
    return s
end

function cleanup_input(varname::Array{String})
    output = String[]
    
    for str in varname
        s = cleanup_input(str)
        push!(output,s)
    end
    return output
end

function isInt(i)
   tol = 1e-3
    return min(i - floor(i), ceil(i)-i) < tol
end

function pick_next(Above, Below)
    (vars_above, vals_above, cost_above, above) = Above
    (vars_below, vals_below, cost_below, below) = Below
    if cost_above == nothing && cost_below == nothing
        println("No integer solution!")
        return fill(nothing, 5)
    elseif cost_above == nothing
        println("Likely infeasible above")
        I_vr, I_vl = identify_INT_vars(vars_below,vals_below)
        return(below, cost_below[1], I_vr, I_vl, "<=")
    elseif cost_below == nothing
        I_vr, I_vl = identify_INT_vars(vars_above,vals_above)
        return(above, cost_above[1], I_vr, I_vl, ">=")
    elseif cost_above[1] <= cost_below[1]
        #println("Chose above")
        I_vr, I_vl = identify_INT_vars(vars_above,vals_above)
        return(above, cost_above[1], I_vr, I_vl, ">=")
    else
        I_vr, I_vl = identify_INT_vars(vars_below,vals_below)
        #println("Chose below")
        return(below, cost_below[1], I_vr, I_vl, "<=")
    end
end

function branch(cur_val,g,best_val, cleaned_vars,direction, pkg_name, mod_name)
    (above, below) = nearest_ints(cur_val)
    Above = g(pkg_name, mod_name,cleaned_vars, [best_val...,above], [direction...,">="])
    Below = g(pkg_name, mod_name,cleaned_vars, [best_val...,below], [direction...,"<="])
    (val, cost, I_vr, I_vl, direct) = pick_next((Above...,above), (Below...,below))

    return (val, cost, I_vr, I_vl, direct)
end

function solve_MIGP(g, pkg_name, model_name, check = true)
    #Solve relaxation
    (vars, vals, cost) = g(pkg_name, model_name)
    GP_count = 1
    I_vr, I_vl = identify_INT_vars(vars,vals)
    I_vr_i = copy(I_vr)
    cleaned_vars = cleanup_input(I_vr)
    
    Nvars = length(I_vr)
    f(x) = 1+(x-1)%Nvars
    direction = fill("", Nvars)
    
    #Initial assignment of all integer variables
    best_cost = Array{Float64}(undef,1,length(I_vr))
    best_val = copy(I_vl)
    for i = 1:length(best_val)
        cur_val = best_val[i]
        if isInt(cur_val)
            best_val[i] = cur_val
            best_cost[i] = cost
            direction[i] = "=="
            deleteat!(I_vr,1)
            deleteat!(I_vl,1)
        else
            (val, cost, I_vr, I_vl, best_dir) = branch(cur_val,g,best_val[1:i-1], cleaned_vars[1:i],direction[1:i-1], pkg_name, model_name)
            direction[i] = best_dir
            best_val[i] = val
            best_cost[i] = cost
            GP_count +=2
        end
    end
    #Repeat to check that there is no change in variables
    j = 1
    repeat = check
    had_change = fill(true, length(best_val))
    max_iter = 100

    while repeat && j < max_iter
        i = f(j)
        cv_copy = copy(cleaned_vars)
        bv_copy = vec(copy(best_val))
        dir_copy = copy(direction)
        #@show i, cv_copy, bv_copy, dir_copy

        #save the old best value
        old_best_value = copy(bv_copy[i])
        var_check = cv_copy[i]
        #Relax one variable
        deleteat!(cv_copy,i)
        deleteat!(bv_copy,i)
        deleteat!(dir_copy,i)

        #@show cv_copy, bv_copy, dir_copy, var_check
        #Solve relaxation of that one variable
        (vars_r, vals_r, cost_r) = g(pkg_name, mod_name,cv_copy, bv_copy, dir_copy)
        GP_count +=1
        
        if vars_r != nothing
            #Get the new optimal value
            I_vr_r, I_vl_r = identify_INT_vars(vars_r,vals_r)
            #@show cleanup_input(I_vr_r), I_vl_r
            #Determine best integer value
            cur_val_i = findall(x -> x == var_check,cleanup_input(I_vr_r))
            cur_val = I_vl_r[cur_val_i[1]]
            #@show cur_val_i, cur_val
            if isInt(cur_val)
                new_best_val = cur_val
                new_dir = "=="
            else
                (new_best_val, new_cost, new_vars, new_vals, new_dir) = branch(cur_val,g,bv_copy, cv_copy,dir_copy, pkg_name, mod_name)
                GP_count +=2
            end
        else
            new_best_val = old_best_value
        end

        if new_best_val != old_best_value
            best_val[i] = new_best_val
            had_change[i] = true
            #@show new_best_val, old_best_value
            #@show new_dir
            direction[i] = new_dir
            #print("Changing "*var_check*"from"*string(old_best_value)*"to "*string(new_best_val))
        else
            had_change[i] = false
        end

        repeat = false
        for h in had_change
            if h
                repeat = true
            end
        end
        j += 1
    end
    return I_vr_i, best_val, cost, GP_count
end
