# Window matheuristic for the single row facility layout problem (SRFLP)
#
# Code for: T. Pammer, M. Sinnl (2026). A matheuristic for solving the single row
# facility layout problem. Computers & Operations Research. https://doi.org/10.1016/j.cor.2026.107397
#
# Usage: see the "Command-line entry point" at the end of this file or the README.

using Random
using DelimitedFiles
using DataFrames
using Gurobi
using JuMP



#read instance with dummy variable for MIP approach
function read_instance_data_dummy_palu(file_path)
    # Read all numbers into a flat array
    data = parse.(Int, split(read(file_path, String)))

    # Step 1: Read n
    n = data[1]

    # Step 2: Read lengths
    l = Vector{Vector{Int64}}()
    for i in 2:n+1
        push!(l, [i-1, data[i]])
    end

    # Step 3: Read condensed upper triangle weights
    num_weights = n*n # number of elements in upper triangle (without diagonal)
    weights_data = data[n+2 : n+1 + num_weights]

    # Step 4: Fill a symmetric matrix
    weights = zeros(Int, n, n)
    idx = 1
    for i in 1:n
        for j in 1:n
            weights[i, j] = weights_data[idx]
            idx += 1
        end
    end

    return n, l, weights
end

quiet(f) = redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        f()
    end
end


#this function uses the 3 dimensional betweennessvariable x from the cplex model and creates the order of the locations in this model
function create_order(x::Array{Float64, 3}, n::Int64)
    
    #create the matrix bij which describes how many elements are between i and j_index
    b = zeros(Float64, n, n)
    for i in 1:n-1
        for j in i+1:n
            b[i,j] = sum(x[i,d,j] for d in 1:n if d != i && d != j)
        end
    end
    i_first, j_last = Tuple(argmax(b))
    b[i_first, i_first] = 1
    for j in 1:i_first-1
        #added two as if there is no object between ifirst and j it is on the second position
        b[i_first, j] = sum(x[j,d,i_first] for d in 1:n if d != i_first && d != j)+2
    end
    for j in i_first+1:n
        b[i_first, j] = sum(x[i_first,d,j] for d in 1:n if d != i_first && d != j)+2
    end
    return b[i_first,:]
    
end





#############################################################################
############################     Main function   ############################
#############################################################################




#function which calculates the costs between all the facilities before and after the window with the facilities in the window for each one seperately
function costs_before_after(cur_sol::Vector{Int64}, window_vars::Vector{Int64}, start_window::Int64, end_window::Int64, costs::Matrix{Int64})
    fac_before = cur_sol[1:start_window-1]
    
    fac_after = cur_sol[end_window+1:length(cur_sol)]
    costs_before = zeros(Int, length(window_vars))
    costs_after = zeros(Int, length(window_vars))
    k=1
    for i in window_vars
        if length(fac_before) >= 1
            costs_before[k] = sum(costs[i,j] for j in fac_before)
        else
            costs_before[k] = 0
        end
        if length(fac_after) >= 1
            costs_after[k] = sum(costs[i,j] for j in fac_after)
        else
            costs_after[k] = 0
        end
        k = k+1
    end
    return costs_before, costs_after
end


function costs_before_window(cur_sol::Vector{Int64}, start_window::Int64, end_window::Int64, costs::Matrix{Int64})
    fac_before = cur_sol[1:start_window-1]
    window_vars = cur_sol[start_window:end_window]
    costs_before = zeros(Int, start_window-1)
    k=1
    for i in fac_before
        if length(fac_before) >= 1
            costs_before[k] = sum(costs[i,j] for j in window_vars)
        else
            costs_before[k] = 0
        end
        k = k+1
    end
    return costs_before
end

function costs_after_window(cur_sol::Vector{Int64}, start_window::Int64, end_window::Int64, costs::Matrix{Int64})
    fac_after = cur_sol[end_window+1:length(cur_sol)]
    window_vars = cur_sol[start_window:end_window]
    costs_after = zeros(Int, length(cur_sol)-end_window)
    k=1
    for i in fac_after
        if length(fac_after) >= 1
            costs_after[k] = sum(costs[i,j] for j in window_vars)
        else
            costs_after[k] = 0
        end
        k = k+1
    end
    return costs_after
end


grb_threads = parse(Int, get(ENV, "GRB_THREADS", "1"))

function MIP_3_dims(len::Vector{Int64}, costs::Matrix{Int64}, size::Int64, fixed_outer_costs::Float64, best_sol_val::Float64, accleration::Bool)
    #create model
    model = quiet(() -> Model(Gurobi.Optimizer))

    # Set the number of threads to 1
    set_optimizer_attribute(model, "Threads", grb_threads)
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "LogToConsole", 0)
    set_optimizer_attribute(model, "LogFile", "")
    #make the variables
    #add the length of 0 for the variables before and after the window
    push!(len, 0)
    push!(len, 0)

    #x[i,k,j] is equal to 1 if location k is between i and j
    @variable(model, x[1:size+2, 1:size+2, 1:size+2], Bin)
    


    K = sum((costs[i,j]*(len[i]+len[j])/2) for i in 1:size+1 for j in i+1:size+2)
    #only until size as facilities n+1 and n+2 are on the outside anyway and so can not be between two other facilities
    @objective(model, Min, K + sum(costs[i,j] * sum(len[k] * x[i,k,j] for k in 1:size if k != j && k != i) for i in 1:size+1 for j in i+1:size+2) + fixed_outer_costs)
    
    #improvement constraint
    if accleration
        @constraint(model, K + sum(costs[i,j] * sum(len[k] * x[i,k,j] for k in 1:size if k != j && k != i) for i in 1:size+1 for j in i+1:size+2) + fixed_outer_costs <= best_sol_val - 0.1)
    end

    @constraint(model, [i=1:size, j = i+1:size+1, k = j+1:size+2],  x[i,j,k] + x[i,k,j] + x[j,i,k] == 1)#7

    @constraint(model, [i=1:size, j = i+1:size+1, k = j+1:size+2, d = 1:size+2; d != i && d != j && d != k], -x[i,d,j] + x[i,d,k] + x[j,d,k] >=0)
    @constraint(model, [i=1:size, j = i+1:size+1, k = j+1:size+2, d = 1:size+2; d != i && d != j && d != k], x[i,d,j] - x[i,d,k] + x[j,d,k] >=0)
    @constraint(model, [i=1:size, j = i+1:size+1, k = j+1:size+2, d = 1:size+2; d != i && d != j && d != k], x[i,d,j] + x[i,d,k] - x[j,d,k] >=0)
    
    @constraint(model, [i=1:size, j = i+1:size+1, k = j+1:size+2, d = 1:size+2; d != i && d != j && d != k], x[i,d,j] + x[i,d,k] + x[j,d,k] <=2)#9

    #add the constraint that the facilities before location and the facilities after location has to be on the outside
    @constraint(model, sum(x[size+1,k,size+2] for k in 1:size) == size)


    quiet(() ->optimize!(model))
    if is_solved_and_feasible(model)
        println(objective_value(model))
        return value.(x)
    else
        return zeros(size,size,size)
    end
    
  
end

#calculate the costs for the facilities out of the window
#in order to check my solutions and to help Christoph
function calc_outer_costs(cur_sol::Vector{Int64}, start_wind::Int64, end_wind::Int64, l::Vector{Int64}, costs::Matrix{Int64}, n::Int64)
    
    #calc the costs for each of the facilities before the window to the facilities in the window
    costs_before_wind = costs_before_window(cur_sol, start_wind, end_wind, costs)
    costs_before = 0
    len_betw = zeros(Int, start_wind-1)
    for i in 1:start_wind-2
        for j in i+1:start_wind-1
            costs_before += costs[cur_sol[i], cur_sol[j]] * ((l[cur_sol[i]]+ l[cur_sol[j]])/ 2 + len_betw[i])
            len_betw[i] += l[cur_sol[j]]
        end
        #add the product of the distance from loc i to the center with the costs between them
        costs_before += costs_before_wind[i] *(len_betw[i]+ l[cur_sol[i]]/2)
    end
    #also add the costs for the facility directly next to the window
    if start_wind > 1
        costs_before += costs_before_wind[start_wind-1] * l[cur_sol[start_wind-1]]/2
    end
    
    len_wind = sum(l[cur_sol[i]] for i in start_wind:end_wind)
    #now calculate the costs between the facilities before and after the window
    for i in 1:start_wind-1
        len_after = 0
        for j in end_wind+1:n
            costs_before += costs[cur_sol[i], cur_sol[j]] * (len_wind+ len_after+ len_betw[i] + (l[cur_sol[i]] + l[cur_sol[j]])/2)
            len_after += l[cur_sol[j]]
        end
    end

    costs_after = 0
    costs_after_wind = costs_after_window(cur_sol, start_wind, end_wind, costs)
    len_betw = l[cur_sol[end_wind+1:n]]
    if end_wind < n
        len_aft = sum(l[cur_sol[i]] for i in end_wind+1:n)
    end
    for i in end_wind+1:n-1
        for j in i+1:n
            costs_after += costs[cur_sol[i], cur_sol[j]] * ((-l[cur_sol[i]]+ l[cur_sol[j]])/ 2 + len_betw[i-end_wind])
            len_betw[i-end_wind] += l[cur_sol[j]]
        end
        #add the product of the distance from loc i to the center with the costs between them
        costs_after += costs_after_wind[i-end_wind] *(len_aft-len_betw[i-end_wind]+ l[cur_sol[i]]/2)
    end
    if end_wind < n
        costs_after += costs_after_wind[n-end_wind] *(len_aft-len_betw[n-end_wind]+ l[cur_sol[n]]/2)
    end
    
    
    return costs_before + costs_after
end

#local search changes 
#using Idea of Palubeckis 2015

#first we need to calculate the so called e values which we need to calc the cut values
function calc_e_values(cur_sol::Vector{Int64}, costs::Matrix{Int64}, n::Int64)
    e_vec = Vector{Int64}(undef, n)
    e_vec[cur_sol[1]] = sum(costs[cur_sol[1], cur_sol[i]] for i in 2:n)
    for m in 2:n-1
        r = cur_sol[m]
        e_vec[r] = sum(costs[r, cur_sol[i]] for i in m+1:n) - sum(costs[r, cur_sol[i]] for i in 1:m-1)
    end
    e_vec[cur_sol[n]] = - sum(costs[cur_sol[n], cur_sol[i]] for i in 1:n-1)
    return e_vec   
end


#calculate the cut values
function calc_cut_values(cur_sol::Vector{Int64}, e_vec::Vector{Int64}, n::Int64)
    cut_vec = Vector{Int64}(undef, n)
    cut_vec[1] = e_vec[cur_sol[1]] # as c0 is zero
    for i in 2:n
        #still have to think about the indexing
        cut_vec[i] = cut_vec[i-1] + e_vec[cur_sol[i]]
    end
   return cut_vec 
end

#calculate cut values when facility r and s gets switched
function calc_cut_prime_values(cur_sol::Vector{Int64}, e_vec::Vector{Int64}, r::Int64, s::Int64, k::Int64, l::Int64, costs::Matrix{Int64}, cut_val_k_1::Int64, cut_vec::Vector{Int64})
    cut_prime_vec = copy(cut_vec[1:l-1])
    cut_prime_vec[k] = cut_val_k_1 + e_vec[s] + 2 * sum(costs[s,cur_sol[i]] for i in k:l-1)
    for m in k+1:l-1
        cut_prime_vec[m] = cut_prime_vec[m-1] + e_vec[cur_sol[m]] + 2*(costs[r,cur_sol[m]] - costs[s,cur_sol[m]])
    end
    return cut_prime_vec
end

function calc_cut_prime_values_insertion(cur_sol::Vector{Int64}, e_vec::Vector{Int64}, r::Int64, k::Int64, l::Int64, costs::Matrix{Int64}, cut_val_k_1::Int64, cut_vec::Vector{Int64})
    if k < l
        cut_prime_vec = copy(cut_vec[1:l-1])
        cut_prime_vec[k] = cut_val_k_1 + e_vec[cur_sol[k+1]] + 2 * costs[r, cur_sol[k+1]]
        for m in k+1:l-1
            cut_prime_vec[m] = cut_prime_vec[m-1] + e_vec[cur_sol[m+1]] + 2 * costs[r, cur_sol[m+1]]
        end
    else
        cut_prime_vec = copy(cut_vec[1:k-1])
        cut_prime_vec[k-1] = cut_val_k_1- e_vec[cur_sol[k-1]] + 2 * costs[r, cur_sol[k-1]]
        for m in 2:k-l
            cut_prime_vec[k-m] = cut_prime_vec[k-m+1] - e_vec[cur_sol[k-m]] + 2*costs[r,cur_sol[k-m]]
        end
    end
    return cut_prime_vec
end



function local_search_insert_phase_palu_new(cur_sol::Vector{Int64}, init_val::Float64, len::Vector{Int64}, costs::Matrix{Int64}, n::Int64, insert_first::Bool, stored_val::Vector{Float64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64}, cut_vec::Vector{Int64})
    gains_matrix = zeros(n,n)
    for k in 1:n
        r = cur_sol[k]
        for l in 1:k-1
            if l > 1
                gains_matrix[k,l] = 2 * b_vec[r,l] + c_minus[k] * (dist[cur_sol[l], cur_sol[k-1]] + (len[cur_sol[l]]+len[cur_sol[k-1]])/2) + len[r]*(cut_vec[l-1]-cut_vec[k-1])
            else
                gains_matrix[k,l] = 2 * b_vec[r,l] + c_minus[k] * (dist[cur_sol[l], cur_sol[k-1]] + (len[cur_sol[l]]+len[cur_sol[k-1]])/2) + len[r]*(0 -cut_vec[k-1])
            end
        end
        for l in k+1:n
            gains_matrix[k,l] = 2 * b_vec[r,l] - c_minus[k] * (dist[cur_sol[k+1], cur_sol[l]] + (len[cur_sol[k+1]]+len[cur_sol[l]])/2) + len[r]*(cut_vec[l]-cut_vec[k])
        end
    end    

    init_sol = copy(cur_sol)
    best_val = init_val
    println("best value before gains calculations: ", best_val)
    println("minimum: ", minimum(gains_matrix))
    while minimum(gains_matrix) < -0.1
        min_index = argmin(gains_matrix)
        # Swap elements to improve the solution
        k = min_index[1]
        l = min_index[2]
        h = cur_sol[k]
        
        
        if k < l
            for i in k:l-1
                cur_sol[i] = cur_sol[i+1]
            end
        else
            for i in 1:k-l 
                cur_sol[k-i+1] = cur_sol[k-i]
            end
        end
        cur_sol[l] = h
        best_val = best_val + gains_matrix[k,l]
        if !insert_first
            local_optima = false
            for val in stored_val
                if val == best_val
                    local_optima = true
                end
            end
            if local_optima
                return init_sol, init_val, false, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
            else
                j1 = min(k,l)
                j2 = max(k,l)
                a_vec = update_a_insert(cur_sol, costs, n, k, l, j1, j2, a_vec)
                b_vec = update_b_insert(cur_sol, len, n, k, l, j1, j2, a_vec, b_vec)
                #maybe change to k and l
                cut_vec = update_c_vec(a_vec, cut_vec, j1,j2,n, cur_sol)
                c_minus = update_c_minus(cut_vec, n, j1, j2, c_minus)
                c_plus = update_c_plus(cut_vec, n, j1, j2, c_plus)
                dist= update_distances(cur_sol, len, n, j1, j2, dist)

                return cur_sol, best_val, true, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
            end
        end
        j1 = min(k,l)
        j2 = max(k,l)
        a_vec = update_a_insert(cur_sol, costs, n, k, l, j1, j2, a_vec)
        b_vec = update_b_insert(cur_sol, len, n, k, l, j1, j2, a_vec, b_vec)
        #maybe change to k and l
        cut_vec = update_c_vec(a_vec, cut_vec, j1,j2,n, cur_sol)
        c_minus = update_c_minus(cut_vec, n, j1, j2, c_minus)
        c_plus = update_c_plus(cut_vec, n, j1, j2, c_plus)
        dist= update_distances(cur_sol, len, n, j1, j2, dist)
        gains_matrix = zeros(n,n)
        for k in 1:n
            r = cur_sol[k]
            for l in 1:k-1
                if l > 1
                    gains_matrix[k,l] = 2 * b_vec[r,l] + c_minus[k] * (dist[cur_sol[l], cur_sol[k-1]] + (len[cur_sol[l]]+len[cur_sol[k-1]])/2) + len[r]*(cut_vec[l-1]-cut_vec[k-1])
                else
                    gains_matrix[k,l] = 2 * b_vec[r,l] + c_minus[k] * (dist[cur_sol[l], cur_sol[k-1]] + (len[cur_sol[l]]+len[cur_sol[k-1]])/2) + len[r]*(0 -cut_vec[k-1])
                end
            end
            for l in k+1:n
                gains_matrix[k,l] = 2 * b_vec[r,l] - c_minus[k] * (dist[cur_sol[k+1], cur_sol[l]] + (len[cur_sol[k+1]]+len[cur_sol[l]])/2) + len[r]*(cut_vec[l]-cut_vec[k])
            end
        end    
    
    end
    if !insert_first
        return cur_sol, best_val, false, a_vec, b_vec, c_minus, c_plus, dist,cut_vec
    end
    local_optima = false
    for val in stored_val
        if val == best_val
            local_optima = true
        end
    end
    if local_optima
        return init_sol, init_val, a_vec, b_vec, c_minus, c_plus, dist,cut_vec
    else 
        return cur_sol, best_val, a_vec, b_vec, c_minus, c_plus, dist,cut_vec
    end
end


function local_search_swap_phase_new(cur_sol::Vector{Int64}, init_val::Float64, len::Vector{Int64}, costs::Matrix{Int64}, n::Int64, insert_first::Bool, stored_val::Vector{Float64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64}, cut_vec::Vector{Int64})
    gains_matrix = zeros(n,n)
    for k in 1:n-1
        for l in k+1:n
            r = cur_sol[k]
            s = cur_sol[l]
            gains_matrix[k,l] = (c_minus[l] - c_minus[k] + 2 * costs[r,s]) * dist[r,s] + 2*(b_vec[r,l-1] + b_vec[s,k+1]) + (len[r] - len[s])/2*(c_plus[l] - c_plus[k])
        end
    end
    
    init_sol = copy(cur_sol)
    best_val = init_val
    i = 0
    while minimum(gains_matrix) < -0.1
        i = i+1
        if i > 100000
            break
        end
        min_index = argmin(gains_matrix)
        k = min_index[1]
        l = min_index[2]
        r = cur_sol[k]
        s = cur_sol[l]
        # Swap elements to improve the solution
        cur_sol[k] = s
        cur_sol[l] = r
        best_val = best_val + gains_matrix[k,l]
        println("best_val after swap: ",best_val)
        r = cur_sol[k]
        s = cur_sol[l]
        if insert_first
            local_optima = false
            for val in stored_val
                if val == best_val
                    local_optima = true
                end
            end
            if local_optima
                return init_sol, init_val, false, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
            else
                a_vec = update_a(cur_sol, costs, n, k, l, a_vec)
                b_vec = update_b(cur_sol, len, n, k, l, a_vec, b_vec)
                #maybe change to k and l
                cut_vec = update_c_vec(a_vec, cut_vec, k, l, n, cur_sol) #update_c_vec(a_vec, zeros(Int64,n), 1,n,n, cur_sol)
                c_minus = update_c_minus(cut_vec, n, k, l, c_minus)
                c_plus = update_c_plus(cut_vec, n, k, l, c_plus)
                dist = update_distances(cur_sol, len, n, k, l, dist)
                return cur_sol, best_val, true, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
            end
        end
        a_vec = update_a(cur_sol, costs, n, k, l, a_vec)
        b_vec = update_b(cur_sol, len, n, k, l, a_vec, b_vec)
        #maybe change to k and l
        cut_vec = update_c_vec(a_vec, cut_vec, k, l, n, cur_sol) #update_c_vec(a_vec, zeros(Int64,n), 1,n,n, cur_sol)
        c_minus = update_c_minus(cut_vec, n, k, l, c_minus)
        c_plus = update_c_plus(cut_vec, n, k, l, c_plus)
        dist = update_distances(cur_sol, len, n, k, l, dist)
        gains_matrix = zeros(n,n)
        for k in 1:n-1
            for l in k+1:n
                r = cur_sol[k]
                s = cur_sol[l]
                gains_matrix[k,l] = (c_minus[l] - c_minus[k] + 2 * costs[r,s]) * dist[r,s] + 2*(b_vec[r,l-1] + b_vec[s,k+1]) + (len[r] - len[s])/2*(c_plus[l] - c_plus[k])
            end
        end
        
    end
    if insert_first
        return cur_sol, best_val, false, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
    end
    local_optima = false
    for val in stored_val
        if val == best_val
            local_optima = true
        end
    end
    if local_optima
        return init_sol, init_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
    else 
        return cur_sol, best_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
    end
end


function mip_improvement_phase(n::Int64, l::Vector{Int64}, costs::Matrix{Int64}, cur_sol::Vector{Int64}, best_sol_val::Float64, best_sol_val_known::Float64, go_back_after_first_impr_MIP::Bool, mult_fac::Float64, stored_val::Vector{Float64}, window_sizes::Vector{Int64}, accleration::Bool)
    if best_sol_val < best_sol_val_known * mult_fac
        MIP_sol = copy(cur_sol)
        value_before_MIP = best_sol_val
        MIP_val = Inf
        loc_impr = false
        h = 1
        while h <= length(window_sizes)
            window_size = window_sizes[h]
            
            for i in 1:floor(Int, n / window_size)
                MIP_sol_before_last_MIP_round = copy(MIP_sol)
                window_vars = MIP_sol[(i - 1) * window_size + 1 : i * window_size]

                # Add the costs before and after to the cost matrix
                costs_before, costs_after = costs_before_after(MIP_sol, window_vars, (i - 1) * window_size + 1, i * window_size, costs)
                costs_window = costs[window_vars, window_vars]
                costs_window = vcat(costs_window, reshape(costs_before, 1, :))
                costs_window = vcat(costs_window, reshape(costs_after, 1, :))
                costs_before = vcat(costs_before, 0, 0)
                costs_after = vcat(costs_after, 0, 0)
                costs_window = hcat(costs_window, reshape(costs_before, :, 1))
                costs_window = hcat(costs_window, reshape(costs_after, :, 1))

                # Fix costs to minimize the objective function
                if accleration
                    min_costs_window_vars_to_outside = zeros(Int, window_size)
                    for j in 1:window_size
                        min_costs_window_vars_to_outside[j] = min(costs_window[window_size + 1, j], costs_window[window_size + 2, j])
                        costs_window[window_size + 1, j] -= min_costs_window_vars_to_outside[j]
                        costs_window[j, window_size + 1] -= min_costs_window_vars_to_outside[j]
                        costs_window[window_size + 2, j] -= min_costs_window_vars_to_outside[j]
                        costs_window[j, window_size + 2] -= min_costs_window_vars_to_outside[j]
                    end

                    # Calculate the fixed costs
                    fixed_costs = calc_outer_costs(MIP_sol, (i - 1) * window_size + 1, i * window_size, l, costs, n)
                    fixed_window_costs = sum(min_costs_window_vars_to_outside) * sum(l[j] for j in window_vars)
                    fixed_costs += fixed_window_costs

                    # Solve the MIP
                    x_MIP = MIP_3_dims(l[window_vars], costs_window, window_size, fixed_costs, best_sol_val, accleration)
                else
                    fixed_costs = calc_outer_costs(MIP_sol, (i - 1) * window_size + 1, i * window_size, l, costs, n)
                    x_MIP = MIP_3_dims(l[window_vars], costs_window, window_size, fixed_costs, best_sol_val, accleration)
                end
                if maximum(x_MIP) > 0
                    println("MIP solved for window ", i, " of size ", window_size)
                    pos_MIP = create_order(x_MIP, window_size)[1:window_size]# .- ones(window_size)
                    help_pos = Vector{Vector{Int64}}()
                    k = 1
                    for loc in window_vars
                        push!(help_pos, [loc, round(Int, pos_MIP[k])])
                        k += 1
                    end

                    # Sort help positions for updating the solution
                    help_pos_rev = copy(help_pos)
                    sort!(help_pos, by = x -> x[2], rev = false)
                    sort!(help_pos_rev, by = x -> x[2], rev = true)
                    help_pos = map(x -> x[1], help_pos)
                    help_pos_rev = map(x -> x[1], help_pos_rev)
                    
                    MIP_sol[(i - 1) * window_size + 1 : i * window_size] = help_pos
                    MIP_sol_rev = copy(MIP_sol)
                    MIP_sol_rev[(i - 1) * window_size + 1 : i * window_size] = help_pos_rev

                    # Calculate the objective value
                    MIP_val = calculate_objective_value_new(MIP_sol, l, n, costs)
                    MIP_val_rev = calculate_objective_value_new(MIP_sol_rev, l, n, costs)

                    if MIP_val_rev < MIP_val
                        MIP_sol = copy(MIP_sol_rev)
                        MIP_val = MIP_val_rev
                    end

                    #add constraint that prevents getting stuck in local optima
                    local_optima = false
                    for val in stored_val
                        if MIP_val == val
                            local_optima = true
                            MIP_val = Inf64
                            MIP_sol = MIP_sol_before_last_MIP_round
                        end
                    end

                    if MIP_val < best_sol_val
                        best_sol_val = MIP_val
                    end

                    #set MIP val to INF so that it is not used



                    if go_back_after_first_impr_MIP
                        if MIP_val < value_before_MIP
                            loc_impr = true
                            best_sol_val = MIP_val
                            cur_sol = MIP_sol
                            value_before_MIP = best_sol_val
                            h = 1
                            return cur_sol, best_sol_val, loc_impr
                            break
                        end
                    end
                end
            end

                

                

            if MIP_val < value_before_MIP
                # If improvement found, restart from the local search phase
                loc_impr = true
                best_sol_val = MIP_val
                cur_sol = MIP_sol
                value_before_MIP = best_sol_val
                h = 1
                return cur_sol, best_sol_val, loc_impr
            else
                h += 1
                loc_impr = false
            end
        end
        return cur_sol, best_sol_val, loc_impr
    else 
        return cur_sol, best_sol_val, false
    end

end


#create a new approach also using scatter search for more diversification
function local_search_MIP_part(n::Int64, l::Vector{Int64}, costs::Matrix{Int64}, best_sol::Vector{Int64}, best_sol_val::Float64, go_back_after_first_impr_MIP::Bool, insert_first::Bool, stored_val::Vector{Float64}, include_MIP::Bool, window_sizes::Vector{Int64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64}, cut_vec::Vector{Int64}, accleration::Bool, time_limit::Int64, start_time::Float64)
    
    loc_impr = true
    
    while loc_impr
        if !insert_first
        #swap phase
            best_sol, best_sol_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec = local_search_swap_phase_new(best_sol, best_sol_val, l, costs, n, insert_first, stored_val, a_vec, b_vec, c_minus, c_plus, dist,cut_vec)

            #insert phase
            best_sol, best_sol_val, improved, a_vec, b_vec, c_minus, c_plus, dist, cut_vec = local_search_insert_phase_palu_new(best_sol, best_sol_val, l[1:n], costs[1:n,1:n], n, insert_first, stored_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec)
            if improved
                continue
            elseif include_MIP
                #mip phase
                #only outcommented because of test of obj value without window approach
                best_sol, best_sol_val, improved = mip_improvement_phase(n, l[1:n], costs[1:n,1:n], best_sol, best_sol_val, best_sol_val, go_back_after_first_impr_MIP, 2.0, stored_val, window_sizes, accleration)
                if improved
                    a_vec = compute_a(best_sol, costs, n)
                    b_vec = compute_b(best_sol, n, a_vec, l)
                    e_vec = calc_e_values(best_sol, costs, n)
                    cut_vec = calc_cut_values(best_sol, e_vec, n)
                    c_minus = compute_c_minus(cut_vec, n)
                    c_plus = compute_c_plus(cut_vec, n)
                    dist = compute_distances(best_sol, l, n)
                    continue
                    
                else
                    loc_impr = false
                end
            else loc_impr = false
            end
        else
            #insert phase
            best_sol, best_sol_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec = local_search_insert_phase_palu_new(best_sol, best_sol_val, l[1:n], costs[1:n,1:n], n, insert_first, stored_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec)
            println("sol after insert: ", best_sol_val)
            #swap phase
            best_sol, best_sol_val, improved, a_vec, b_vec, c_minus, c_plus, dist, cut_vec = local_search_swap_phase_new(best_sol, best_sol_val, l, costs, n, insert_first, stored_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec)
            
            if improved
                continue
            elseif include_MIP
                println("solution after ls phase: ", best_sol_val)
                if time() - start_time > time_limit - (n/10)
                    return best_sol, best_sol_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
                end
                #mip phase
                best_sol, best_sol_val, improved = mip_improvement_phase(n, l[1:n], costs[1:n,1:n], best_sol, best_sol_val, best_sol_val, go_back_after_first_impr_MIP, 2.0, stored_val, window_sizes, accleration)
                if improved
                    println("improved by MIP!")
                    a_vec = compute_a(best_sol, costs, n)
                    b_vec = compute_b(best_sol, n, a_vec, l)
                    e_vec = calc_e_values(best_sol, costs, n)
                    cut_vec = calc_cut_values(best_sol, e_vec, n)
                    c_minus = compute_c_minus(cut_vec, n)
                    c_plus = compute_c_plus(cut_vec, n)
                    dist = compute_distances(best_sol, l, n)
                    continue
                else
                    loc_impr = false
                end
            else
                loc_impr = false
            end
        end
    end
    return best_sol, best_sol_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec
end

function compute_T_max(cur_sol::Vector{Int64}, n::Int64, costs::Matrix{Int64}, len::Vector{Int64}, sample_size::Int64)
    #generate necessary data    
    a_vec = compute_a(cur_sol, costs, n)
    b_vec = compute_b(cur_sol, n, a_vec, len)
    cut_vec = update_c_vec(a_vec, zeros(Int64,n), 1, n, n, cur_sol)
    c_minus = compute_c_minus(cut_vec, n)
    c_plus = compute_c_plus(cut_vec, n)
    dist = compute_distances(cur_sol, len, n)
    T_max = 0
    i = 1
    while i < sample_size
        k = rand(1:n-1)
        l = rand(k+1:n)
        r = cur_sol[k]
        s = cur_sol[l]
        move_gain = (c_minus[l] - c_minus[k] + 2 * costs[r,s]) * dist[r,s] + 2*(b_vec[r,l-1] + b_vec[s,k+1]) + (len[r] - len[s])/2*(c_plus[l] - c_plus[k])
        if move_gain < 0
            move_gain = move_gain * -1
        end
        if move_gain > T_max
            T_max = move_gain
        end
        i += 1
    end
    return T_max
    
end

function move_gain_l_k_1(e_vec::Vector{Int64}, cut_vec::Vector{Int64}, costs::Matrix{Int64}, len::Vector{Int64}, r::Int64, k::Int64, s::Int64, l::Int64)
    move_gain = 0
    cut_prime_value_k = 0
    if k > 1
        move_gain = (cut_vec[k-1] - cut_vec[l]) * (len[s]-len[r])/2 + (cut_vec[k-1] + e_vec[s]+ 2*costs[r,s] - cut_vec[k]) * (len[s]+len[r])/2
        cut_prime_value_k = cut_vec[k-1]+ e_vec[s] + 2 * costs[s,r]
    else
        move_gain = (-cut_vec[l]) * (len[s]-len[r])/2 + (e_vec[s]+ 2*costs[r,s] - cut_vec[k]) * (len[s]+len[r])/2
        cut_prime_value_k = 0+ e_vec[s] + 2 * costs[s,r]
    end 
    return move_gain, cut_prime_value_k   
end

function compute_distances(cur_sol::Vector{Int64}, len::Vector{Int64}, n::Int64)
    d = zeros(Float64, n,n)
    for i in 1:n-1
        k = cur_sol[i]
        dist_help = len[k]/2
        for j in i+1:n
            l = cur_sol[j]
            d[k,l] = dist_help + len[l]/2
            d[l,k] = d[k,l]
            dist_help += len[l]
        end
    end
    return d
end

function calculate_objective_value_new(cur_sol::Vector{Int64}, len::Vector{Int64}, n::Int64, costs::Matrix{Int64})
    wd = zeros(Float64, n,n)
    sum_wd = 0.0
    for i in 1:n-1
        k = cur_sol[i]
        dist_help = len[k]/2
        for j in i+1:n
            l = cur_sol[j]
            wd[k,l] = (dist_help + len[l]/2) * costs[k,l]
            sum_wd += wd[k,l]
            dist_help += len[l]
        end
    end
    return sum_wd    
end

function update_distances(cur_sol::Vector{Int64}, len::Vector{Int64}, n::Int64, k::Int64, l::Int64, dist::Matrix{Float64})
    for m in 1:k-1
        u = cur_sol[m]
        for j = k:l
            dist[u,cur_sol[j]] = dist[u,cur_sol[j-1]] + (len[cur_sol[j]] + len[cur_sol[j-1]])/2
            dist[cur_sol[j],u] = dist[u,cur_sol[j]]
        end
    end
    for m in k:l
        u = cur_sol[m]
        if m == k || m == l
            j1 = m-1
            j2 = m+1
        else
            j1 = k
            j2 = l
        end
        for j in j2:n
            dist[u,cur_sol[j]] = dist[u,cur_sol[j-1]] + (len[cur_sol[j]] + len[cur_sol[j-1]])/2
            dist[cur_sol[j],u] = dist[u,cur_sol[j]]
        end
    end
    return dist
    
end

function compute_a(cur_sol::Vector{Int64}, costs::Matrix{Int64}, n::Int64)
    a = zeros(Int64, n,n)
    for m in 1:n
        u = cur_sol[m]
        for j in m+1:n
            a[u,j] = a[u,j-1] + costs[u, cur_sol[j]] 
        end
        for j in 1:m-1
            a[u,m-j] = a[u,m-j+1] + costs[u, cur_sol[m-j]]
        end 
    end
    return a
end

function update_a(cur_sol::Vector{Int64}, costs::Matrix{Int64}, n::Int64, k::Int64, l::Int64, a_vec::Matrix{Int64})
    for m in 1:k-1
        u = cur_sol[m]
        for j = k:l-1
            a_vec[u,j] = a_vec[u,j-1] + costs[u, cur_sol[j]]
        end
    end
    for m in l+1:n
        u = cur_sol[m]
        for j in 1:l-k 
            a_vec[u,l-j+1] = a_vec[u,l-j+2] + costs[u, cur_sol[l-j+1]]
        end
    end
    for m in k:l
        u = cur_sol[m]
        if m == k || m == l
            j1 = m-1
            j2 = m+1
        else
            j1 = k
            j2 = l
        end
        a_vec[u,m] = 0
        for j in 1:j1
            a_vec[u,j1-j+1] = a_vec[u,j1-j+2] + costs[u, cur_sol[j1-j+1]]
        end
        for j in j2:n 
            a_vec[u,j] = a_vec[u,j-1] + costs[u, cur_sol[j]]
        end
    end
    return a_vec    
end

function update_a_insert(cur_sol::Vector{Int64}, costs::Matrix{Int64}, n::Int64, k::Int64, l::Int64, j1::Int64, j2::Int64, a_vec::Matrix{Int64})
    for m in 1:j1-1
        u = cur_sol[m]
        for j = j1:j2-1
            a_vec[u,j] = a_vec[u,j-1] + costs[u, cur_sol[j]]
        end
    end
    for m in j2+1:n
        u = cur_sol[m]
        for j in 1:j2-j1
            a_vec[u,j2-j+1] = a_vec[u,j2-j+2] + costs[u, cur_sol[j2-j+1]]
        end
    end
    
    #4. m == l
    u = cur_sol[l]
    a_vec[u,l] = 0
    for j in 1:l-1
        a_vec[u, l-j] = a_vec[u,l-j+1] + costs[u, cur_sol[l-j]]
    end
    for j in l+1:n
        a_vec[u, j] = a_vec[u,j-1] + costs[u, cur_sol[j]]
    end
    if k < l
        for m in j1:j2-1
            u = cur_sol[m]
            for j in j1:j2-1
                a_vec[u,j] = a_vec[u,j+1]
            end
            for j in 1:j1-1
                a_vec[u,j1-j] = a_vec[u, j1-j+1] + costs[u, cur_sol[j1-j]]
            end
            for j in j2:n
                a_vec[u, j] = a_vec[u,j-1] + costs[u, cur_sol[j]]
            end
        end
    else
        for m in j1+1:j2
            u = cur_sol[m]
            for j in 1:j2-j1
                a_vec[u,j2-j+1] = a_vec[u, j2-j] 
            end
            for j in 1:j1
                a_vec[u,j1-j+1] = a_vec[u, j1-j+2] + costs[u, cur_sol[j1-j+1]]
            end
            for j in j2+1:n
                a_vec[u, j] = a_vec[u,j-1] + costs[u, cur_sol[j]]
            end
        end
    end
    return a_vec    
end

function compute_b(cur_sol::Vector{Int64}, n::Int64, a::Matrix{Int64}, len::Vector{Int64})
    b = zeros(Float64, n,n)
    for m in 1:n
        helper = 0
        u = cur_sol[m]
        for j in m+1:n
            if j > m+1
                helper += a[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
            end
            b[u,j] = a[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
        end
        helper = 0
        for j in 1:m-1
            if j > 1
                helper += a[u,m-j+1] * (len[cur_sol[m-j]] + len[cur_sol[m-j+1]])/2
            end
            b[u,m-j] = a[u,m-j] * (len[u] + len[cur_sol[m-j]])/2  + helper
            
        end
    end
    return b   
end

function update_b(cur_sol::Vector{Int64}, len::Vector{Int64}, n::Int64, k::Int64, l::Int64, a_vec::Matrix{Int64}, b_vec::Matrix{Float64})
    for m in 1:k-1
        u = cur_sol[m]
        helper = 0
        for j in m+2:k-1
            helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
        end
        for j = k:n
            if j > m+1
                helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
            end
            b_vec[u,j] = a_vec[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
        end
    end

    for m in l+1:n
        u = cur_sol[m]
        helper = 0
        for j in 1:m-l
            if j > 1
                helper += a_vec[u,m-j+1] * (len[cur_sol[m-j]] + len[cur_sol[m-j+1]])/2
            end
        end    
        for j in 1:l
            if j > 1
                helper += a_vec[u,l-j+2] * (len[cur_sol[l-j+1]] + len[cur_sol[l-j+2]])/2
            end
            b_vec[u,l-j+1] = a_vec[u,l-j+1] * (len[u] + len[cur_sol[l-j+1]])/2  + helper
        end
    end
    for m in k:l
        u = cur_sol[m]
        if m == k || m == l
            j1 = m-1
            j2 = m+1
        else
            j1 = k
            j2 = l
        end
        b_vec[u,m] = 0
        helper = 0
        for j in 1:m-j1
            if j > 1
                helper += a_vec[u,m-j+1] * (len[cur_sol[m-j]] + len[cur_sol[m-j+1]])/2
            end
        end   
        for j in 1:j1
            if j > 1
                helper += a_vec[u,j1-j+2] * (len[cur_sol[j1-j+1]] + len[cur_sol[j1-j+2]])/2
            end
            b_vec[u,j1-j+1] = a_vec[u,j1-j+1] * (len[u] + len[cur_sol[j1-j+1]])/2  + helper
        end
        helper = 0
        for j in m+2:j2-1
            helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
        end
        for j in j2:n
            if j > m+1
                helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
            end
            b_vec[u,j] = a_vec[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
        end
    end
    return b_vec    
end

function update_b_insert(cur_sol::Vector{Int64}, len::Vector{Int64}, n::Int64, k::Int64, l::Int64, j1::Int64, j2::Int64, a_vec::Matrix{Int64}, b_vec::Matrix{Float64})
    for m in 1:j1-1
        u = cur_sol[m]
        helper = 0
        for j in m+2:j1-1
            helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
        end
        for j = j1:n
            if j > m+1
                helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
            end
            b_vec[u,j] = a_vec[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
        end
    end

    for m in j2+1:n
        u = cur_sol[m]
        helper = 0
        for j in 1:m-j2
            if j > 1
                helper += a_vec[u,m-j+1] * (len[cur_sol[m-j]] + len[cur_sol[m-j+1]])/2
            end
        end    
        for j in 1:j2
            if j > 1
                helper += a_vec[u,j2-j+2] * (len[cur_sol[j2-j+1]] + len[cur_sol[j2-j+2]])/2
            end
            b_vec[u,j2-j+1] = a_vec[u,j2-j+1] * (len[u] + len[cur_sol[j2-j+1]])/2  + helper
        end
    end

    #4. m == l
    u = cur_sol[l]
    b_vec[u,l] = 0
    helper = 0
    for j in 1:l-1
        if j > 1
            helper += a_vec[u,l-j+1] * (len[cur_sol[l-j]] + len[cur_sol[l-j+1]])/2
        end
        b_vec[u,l-j] = a_vec[u,l-j] * (len[u] + len[cur_sol[l-j]])/2  + helper
    end
    helper = 0
    for j in l+1:n
        if j > l+1
            helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
        end
        b_vec[u,j] = a_vec[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
    end

    if k < l
        #we can replace all k with j1 and outcomment this if condition
        #we only also need to set the current for loops with j1 to j1-1
        if k > 1
            j1 = k-1
        else
            j1 = k 
        end
        j2 = l
        for m in k:l-1
            u = cur_sol[m]
            for j in k:l-1
                b_vec[u,j] = b_vec[u,j+1]
            end
            helper = 0
            for j in m+2:j2-1
                helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
            end
            for j = j2:n
                if j > m+1
                    helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
                end
                b_vec[u,j] = a_vec[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
            end
            helper = 0
            for j in 2:m-j1
                helper += a_vec[u,m-j+1] * (len[cur_sol[m-j]] + len[cur_sol[m-j+1]])/2
            end
            for j in 1:j1
                if j > 1
                    helper += a_vec[u,j1-j+2] * (len[cur_sol[j1-j+1]] + len[cur_sol[j1-j+2]])/2
                end
                b_vec[u,j1-j+1] = a_vec[u,j1-j+1] * (len[u] + len[cur_sol[j1-j+1]])/2  + helper
            end
        end
    else
        j1 = l
        if k < n
            j2 = k+1
        else
            j2 = k
        end
        for m in l+1:k
            u = cur_sol[m]
            for j in 1:k-l
                b_vec[u,k-j+1] = b_vec[u, k-j] 
            end
            helper = 0
            for j in m+2:j2-1
                helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
            end
            for j = j2:n
                if j > m+1
                    helper += a_vec[u,j-1] * (len[cur_sol[j-1]] + len[cur_sol[j]])/2
                end
                b_vec[u,j] = a_vec[u,j] * (len[u] + len[cur_sol[j]])/2 + helper
            end
            helper = 0
            for j in 2:m-j1
                helper += a_vec[u,m-j+1] * (len[cur_sol[m-j]] + len[cur_sol[m-j+1]])/2
            end
            for j in 1:j1
                if j > 1
                    helper += a_vec[u,j1-j+2] * (len[cur_sol[j1-j+1]] + len[cur_sol[j1-j+2]])/2
                end
                b_vec[u,j1-j+1] = a_vec[u,j1-j+1] * (len[u] + len[cur_sol[j1-j+1]])/2  + helper
            end
        end
    end
    
    return b_vec    
end

function compute_c_minus(cut_vec::Vector{Int64}, n::Int64)
    c_minus = zeros(Int64,n)
    c_minus[1] = cut_vec[1]
    for i in 2:n-1
        c_minus[i] = cut_vec[i]-cut_vec[i-1]
    end
    c_minus[n] = -cut_vec[n-1]
    return c_minus
end

function update_c_minus(cut_vec::Vector{Int64}, n::Int64, k::Int64, l::Int64, c_minus::Vector{Int64})
    if k == 1
        c_minus[1] = cut_vec[1]
        k = 2
    end
    if l == n
        c_minus[n] = -cut_vec[n-1]
        l = n-1
    end
    for i in k:l+1
        c_minus[i] = cut_vec[i] - cut_vec[i-1]
    end
    return c_minus
end

function compute_c_plus(cut_vec::Vector{Int64}, n::Int64)
    c_plus = zeros(Int64,n)
    c_plus[1] = cut_vec[1]
    for i in 2:n-1
        c_plus[i] = cut_vec[i]+cut_vec[i-1]
    end
    c_plus[n] = cut_vec[n-1]
    return c_plus
end

function update_c_plus(cut_vec::Vector{Int64}, n::Int64, k::Int64, l::Int64, c_plus::Vector{Int64})
    if k == 1
        c_plus[1] = cut_vec[1]
        k = 2
    end
    if l == n
        c_plus[n] = cut_vec[n-1]
        l = n-1
    end
    for i in k:l+1
        c_plus[i] = cut_vec[i] + cut_vec[i-1]
    end
    return c_plus
end

function update_c_vec(a::Matrix{Int64}, cut_vec::Vector{Int64}, k::Int64, l::Int64, n::Int64, cur_sol::Vector{Int64})
    if k == 1
        cut_vec[1] = a[cur_sol[1], n]
        for i in k+1:l-1
            cut_vec[i] = cut_vec[i-1] + a[cur_sol[i], n] - a[cur_sol[i], 1]
        end
    else
        for i in k:l-1
            cut_vec[i] = cut_vec[i-1] + a[cur_sol[i], n] - a[cur_sol[i], 1]
        end
    end
    return cut_vec    
end

function generate_interchange(cur_sol::Vector{Int64}, beta::Int64, beta_prime::Int64, z::Int64, gain_com_prob::Float64, r::Int64, s::Int64, k::Int64, l::Int64, e_vec::Vector{Int64}, cut_vec::Vector{Int64}, n::Int64, costs::Matrix{Int64}, len::Vector{Int64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64})
    move_gain = 0
    cut_prime_values = Vector{Int64}()
    if beta <= round.(Int64, gain_com_prob * beta_prime)
        if beta == 1 && z == 1
            e_vec = calc_e_values(cur_sol, costs, n)
            cut_vec = calc_cut_values(cur_sol, e_vec, n)
        end
        if l > k+1
            if k >1
                cut_prime_values = calc_cut_prime_values(cur_sol, e_vec, r, s, k, l, costs, cut_vec[k-1], cut_vec)
                move_gain = (cut_vec[k-1]- cut_vec[l])*((len[s]-len[r])/2) + cut_prime_values[k]* (len[s]+len[cur_sol[k+1]])/2- cut_vec[k]*(len[r] + len[cur_sol[k+1]])/2+ cut_prime_values[l-1] *(len[cur_sol[l-1]]+ len[r])/2- cut_vec[l-1] * (len[cur_sol[l-1]] + len[s])/2
                if l > k+2
                    move_gain += sum((cut_prime_values[m] - cut_vec[m]) * (len[cur_sol[m]] + len[cur_sol[m+1]])/2 for m in k+1:l-2)
                end
            else
                cut_prime_values = calc_cut_prime_values(cur_sol, e_vec, r, s, k, l, costs, 0, cut_vec)
            
                move_gain =  -cut_vec[l]*(len[s]-len[r])/2 + cut_prime_values[k]* (len[s]+len[cur_sol[k+1]])/ 2 - cut_vec[k]*(len[r] + len[cur_sol[k+1]])/2 + cut_prime_values[l-1] *(len[cur_sol[l-1]]+ len[r])/2 - cut_vec[l-1] * (len[cur_sol[l-1]]+ len[s])/2
                if l > k+2
                    move_gain += sum((cut_prime_values[m] - cut_vec[m]) * (len[cur_sol[m]] + len[cur_sol[m+1]])/2 for m in k+1:l-2)
                end
            end


        else
            cut_prime_values = copy(cut_vec[1:l-1])
            move_gain, cut_prime_values[k] = move_gain_l_k_1(e_vec, cut_vec, costs, len, r, k, s, l)
        end
    else
        if beta == round.(Int64, gain_com_prob * beta_prime) + 1 && z == 1
            a_vec = compute_a(cur_sol, costs, n)
            b_vec = compute_b(cur_sol, n, a_vec, len)
            if round.(Int64, gain_com_prob * beta_prime) == 0
                cut_vec = update_c_vec(a_vec, zeros(Int64,n), 1, n, n, cur_sol)
            end
            c_minus = compute_c_minus(cut_vec, n)
            c_plus = compute_c_plus(cut_vec, n)
            dist = compute_distances(cur_sol, len, n)
        end
        move_gain = (c_minus[l] - c_minus[k] + 2 * costs[r,s]) * dist[r,s] + 2*(b_vec[r,l-1] + b_vec[s,k+1]) + (len[r] - len[s])/2*(c_plus[l] - c_plus[k])

    end
    return move_gain, e_vec, cut_vec, cut_prime_values, a_vec, b_vec, c_minus, c_plus, dist
end

function perform_interchange(cur_sol::Vector{Int64}, r::Int64, k::Int64, s::Int64, l::Int64, beta::Int64, beta_prime::Int64, gain_com_prob::Float64, e_vec::Vector{Int64}, cut_vec::Vector{Int64}, cut_prime_vec::Vector{Int64}, n::Int64, costs::Matrix{Int64}, len::Vector{Int64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64})
    #swap positions in the current solution
    cur_sol[k] = s
    cur_sol[l] = r
    if beta <= round.(Int64, gain_com_prob * beta_prime)
        #update the vectors
        cut_vec[k:l-1] = cut_prime_vec[k:l-1]
        #we multiply by two as it was reduced before by that amount
        for i in k:l-1
            if i > 1
                e_vec[cur_sol[i]] = cut_prime_vec[i] - cut_prime_vec[i-1]
            else
                e_vec[cur_sol[i]] = cut_prime_vec[i]
            end
        end
        e_vec[r] = e_vec[r] - 2 * sum(costs[r, cur_sol[i]] for i in k:l-1)
    #still need to find out how the else condition works
    else
        a_vec = update_a(cur_sol, costs, n, k, l, a_vec)
        b_vec = update_b(cur_sol, len, n, k, l, a_vec, b_vec)
        #maybe change to k and l
        cut_vec = update_c_vec(a_vec, cut_vec, k, l, n, cur_sol) #update_c_vec(a_vec, zeros(Int64,n), 1,n,n, cur_sol)
        c_minus = update_c_minus(cut_vec, n, k, l, c_minus)
        c_plus = update_c_plus(cut_vec, n, k, l, c_plus)
        dist = update_distances(cur_sol, len, n, k, l, dist)
    end
    return cur_sol, e_vec, cut_vec, a_vec, b_vec, c_minus, c_plus, dist
end


function generate_insertion(cur_sol::Vector{Int64}, beta::Int64, beta_prime::Int64, z::Int64, gain_com_prob::Float64, r::Int64, k::Int64, l::Int64, e_vec::Vector{Int64}, cut_vec::Vector{Int64}, n::Int64, costs::Matrix{Int64}, len::Vector{Int64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64})
    move_gain = 0
    cut_prime_values = Vector{Int64}()
    if beta <= round.(Int64, gain_com_prob * beta_prime)
        if beta == 1 && z == 1
            e_vec = calc_e_values(cur_sol, costs, n)
            cut_vec = calc_cut_values(cur_sol, e_vec, n)
        end
        if k < l
            if l > k+1
                if k >1
                    cut_prime_values = calc_cut_prime_values_insertion(cur_sol, e_vec, r, k, l, costs, cut_vec[k-1], cut_vec)
                    move_gain = cut_vec[k-1]*((len[cur_sol[k+1]]-len[r])/2) - cut_vec[k]*(len[r] + len[cur_sol[k+1]])/2+ cut_prime_values[l-1] *(len[cur_sol[l]]+ len[r])/2 + cut_vec[l] * (-len[cur_sol[l]] + len[r])/2
                    move_gain += sum((cut_prime_values[m-1] - cut_vec[m]) * (len[cur_sol[m]] + len[cur_sol[m+1]])/2 for m in k+1:l-1)
                else
                    cut_prime_values = calc_cut_prime_values_insertion(cur_sol, e_vec, r, k, l, costs, 0, cut_vec)
                    move_gain = - cut_vec[k]*(len[r] + len[cur_sol[k+1]])/2+ cut_prime_values[l-1] *(len[cur_sol[l]]+ len[r])/2 + cut_vec[l] * (-len[cur_sol[l]] + len[r])/2
                    move_gain += sum((cut_prime_values[m-1] - cut_vec[m]) * (len[cur_sol[m]] + len[cur_sol[m+1]])/2 for m in k+1:l-1)
                end


            else
                cut_prime_values = copy(cut_vec[1:l-1])
                s = cur_sol[l]
                move_gain, cut_prime_values[k] = move_gain_l_k_1(e_vec, cut_vec, costs, len, r, k, s, l)
            end
        else
            if k > l+1
                if l >1
                    cut_prime_values = calc_cut_prime_values_insertion(cur_sol, e_vec, r, k, l, costs, cut_vec[k], cut_vec)
                    move_gain = cut_vec[k]*((len[cur_sol[k-1]]-len[r])/2) - cut_vec[k-1]*(len[r] + len[cur_sol[k-1]])/2+ cut_prime_values[l] *(len[cur_sol[l]]+ len[r])/2 + cut_vec[l-1] * (-len[cur_sol[l]] + len[r])/2
                    move_gain += sum((cut_prime_values[m+1] - cut_vec[m]) * (len[cur_sol[m]] + len[cur_sol[m+1]])/2 for m in l:k-2)
                else
                    cut_prime_values = calc_cut_prime_values_insertion(cur_sol, e_vec, r, k, l, costs, cut_vec[k], cut_vec)
                    move_gain = cut_vec[k]*((len[cur_sol[k-1]]-len[r])/2) - cut_vec[k-1]*(len[r] + len[cur_sol[k-1]])/2+ cut_prime_values[l] *(len[cur_sol[l]]+ len[r])/2
                    move_gain += sum((cut_prime_values[m+1] - cut_vec[m]) * (len[cur_sol[m]] + len[cur_sol[m+1]])/2 for m in l:k-2)
                end


            else
                cut_prime_values = copy(cut_vec[1:k-1])
                s = cur_sol[l]
                #parameters are switched here as l smaller k
                move_gain, cut_prime_values[l] = move_gain_l_k_1(e_vec, cut_vec, costs, len, s, l, r, k)
            end
        end
    else
        #implement later on
        if beta == round.(Int64, gain_com_prob * beta_prime) + 1 && z == 1
            a_vec = compute_a(cur_sol, costs, n)
            b_vec = compute_b(cur_sol, n, a_vec, len)
            if round.(Int64, gain_com_prob * beta_prime) == 0
                cut_vec = update_c_vec(a_vec, zeros(Int64,n), 1, n, n, cur_sol)
            end
            c_minus = compute_c_minus(cut_vec, n)
            c_plus = compute_c_plus(cut_vec, n)
            dist = compute_distances(cur_sol, len, n)
        end
        if k < l
            move_gain = 2 * b_vec[r,l] - c_minus[k] * (dist[cur_sol[k+1], cur_sol[l]] + (len[cur_sol[k+1]]+len[cur_sol[l]])/2) + len[r]*(cut_vec[l]-cut_vec[k])
        else
            if l > 1
                move_gain = 2 * b_vec[r,l] + c_minus[k] * (dist[cur_sol[l], cur_sol[k-1]] + (len[cur_sol[l]]+len[cur_sol[k-1]])/2) + len[r]*(cut_vec[l-1]-cut_vec[k-1])
            else
                move_gain = 2 * b_vec[r,l] + c_minus[k] * (dist[cur_sol[l], cur_sol[k-1]] + (len[cur_sol[l]]+len[cur_sol[k-1]])/2) + len[r]*(0 -cut_vec[k-1])
            end
        end    
    end
    return move_gain, e_vec, cut_vec, cut_prime_values, a_vec, b_vec, c_minus, c_plus, dist
end

function perform_insertion(cur_sol::Vector{Int64}, r::Int64, k::Int64, l::Int64, beta::Int64, beta_prime::Int64, gain_com_prob::Float64, e_vec::Vector{Int64}, cut_vec::Vector{Int64}, cut_prime_vec::Vector{Int64}, n::Int64, costs::Matrix{Int64}, len::Vector{Int64}, a_vec::Matrix{Int64}, b_vec::Matrix{Float64}, c_minus::Vector{Int64}, c_plus::Vector{Int64}, dist::Matrix{Float64})
    if beta <= round.(Int64, gain_com_prob * beta_prime)
        if k <l
            #update positions
            for i in k:l-1
                cur_sol[i] = cur_sol[i+1]
            end
            cur_sol[l] = r
            #update the vectors
            cut_vec[k:l-1] = cut_prime_vec[k:l-1]
            #we multiply by two as it was reduced before by that amount
            for i in k:l-1
                #can be calculated that way because cur_sol already updated
                e_vec[cur_sol[i]] += 2* costs[r,cur_sol[i]] 
            end
            e_vec[r] -= 2 * sum(costs[r, cur_sol[i]] for i in k:l-1)
        else
            #update positions
            for i in 1:k-l
                cur_sol[k-i+1] = cur_sol[k-i]
            end
            cur_sol[l] = r
            #update the vectors
            cut_vec[l:k-1] = cut_prime_vec[l:k-1]
            #we multiply by two as it was reduced before by that amount
            for i in l+1:k
                e_vec[cur_sol[i]] -= 2* costs[r,cur_sol[i]]
            end
            e_vec[r] += 2 * sum(costs[r, cur_sol[i]] for i in l+1:k)
        end
    #still need to find out how the else condition works
    else
        if k <l
            #update positions
            for i in k:l-1
                cur_sol[i] = cur_sol[i+1]
            end
            cur_sol[l] = r
        else
            #update positions
            for i in 1:k-l
                cur_sol[k-i+1] = cur_sol[k-i]
            end
            cur_sol[l] = r
        end
        j1 = min(k,l)
        j2 = max(k,l)
        a_vec = update_a_insert(cur_sol, costs, n, k, l, j1, j2, a_vec)
        b_vec = update_b_insert(cur_sol, len, n, k, l, j1, j2, a_vec, b_vec)
        #maybe change to k and l
        cut_vec = update_c_vec(a_vec, cut_vec, j1,j2,n, cur_sol)
        c_minus = update_c_minus(cut_vec, n, j1, j2, c_minus)
        c_plus = update_c_plus(cut_vec, n, j1, j2, c_plus)
        dist= update_distances(cur_sol, len, n, j1, j2, dist)
    end
    return cur_sol, e_vec, cut_vec, a_vec, b_vec, c_minus, c_plus, dist

end

function MSA(cooling_fac::Float64, n_inner_loop::Int64, T_min::Float64, gain_com_prob::Float64, prob_interchange::Float64, n::Int64, costs::Matrix{Int64}, len::Vector{Vector{Int64}}, accleration::Bool, time_limit::Int64)
    len = [leng[2] for leng in len]
    cur_sol = randperm(n)
    cur_sol_val = calculate_objective_value_new(cur_sol, len, n, costs)
    best_sol = copy(cur_sol)
    best_sol_val = cur_sol_val
    best_sol_bef_MIP_val = cur_sol_val
    best_ovr_sol_val = best_sol_val
    best_ovr_sol = copy(best_sol)
    T_max = compute_T_max(best_sol, n, costs, len, 5000)
    beta_prime = round.(Int64, (log(T_min) - log(T_max))/log(cooling_fac))
    println("Number outer loops: ", beta_prime)
    n_wo_impr = 0
    accept = true
    e_vec = Vector{Int64}()
    cut_vec = Vector{Int64}()
    cut_prime_vec = Vector{Int64}()
    a_vec = Matrix{Int64}(undef, 0,0)
    b_vec = Matrix{Float64}(undef, 0,0)
    c_minus = Vector{Int64}()
    c_plus = Vector{Int64}()
    dist = Matrix{Float64}(undef, 0,0)
    #including time constraint for the loop
    start_time = time()
    n_SA_starts = 0
    while (time() - start_time) < time_limit  - (n/2)
        n_SA_starts += 1
        best_sol_val = cur_sol_val
        T = T_max
        for beta in 1:beta_prime
            for z in 1:n_inner_loop
                rand_num = rand()
                k = 0
                r = 0
                l = 0
                s = 0
                if rand_num <= prob_interchange
                    #randomly choose two locations k and l where k < l with their corresponding facilities r and s
                    k = rand(1:n-1)
                    r = cur_sol[k]
                    l = rand(k+1:n)
                    s = cur_sol[l]
                    move_gain, e_vec, cut_vec, cut_prime_vec, a_vec, b_vec, c_minus, c_plus, dist = generate_interchange(cur_sol, beta, beta_prime, z, gain_com_prob, r, s, k, l, e_vec, cut_vec, n, costs, len, a_vec, b_vec, c_minus, c_plus, dist)
                else
                    #randomly choose two locations k and l with k's corresponding facilities r 
                    k = rand(1:n)
                    r = cur_sol[k]
                    l = rand(1:n)
                    while k == l 
                        l = rand(1:n)
                    end
                    move_gain, e_vec, cut_vec, cut_prime_vec, a_vec, b_vec, c_minus, c_plus, dist = generate_insertion(cur_sol, beta, beta_prime, z, gain_com_prob, r, k, l, e_vec, cut_vec, n, costs, len, a_vec, b_vec, c_minus, c_plus, dist)
                end
                if move_gain <= 0
                    accept = true
                else
                    acc_prob = rand()
                    #here the acceptance of a new solution is dependent on the coolin factor and the gain value
                    #the higher the temperature the higher the chance for diversification
                    if acc_prob <= exp(-move_gain/T)
                        accept = true
                    else
                        accept = false
                    end
                end
                if accept
                    cur_sol_val += move_gain
                    if rand_num <= prob_interchange
                        cur_sol, e_vec, cut_vec, a_vec, b_vec, c_minus, c_plus, dist = perform_interchange(cur_sol, r,k , s, l, beta, beta_prime, gain_com_prob, e_vec, cut_vec, cut_prime_vec, n, costs, len, a_vec, b_vec, c_minus, c_plus, dist)
                    else
                        cur_sol, e_vec, cut_vec, a_vec, b_vec, c_minus, c_plus, dist = perform_insertion(cur_sol, r, k, l, beta, beta_prime, gain_com_prob, e_vec, cut_vec, cut_prime_vec, n, costs, len, a_vec, b_vec, c_minus, c_plus, dist)
                    end
                    if cur_sol_val < best_sol_val
                        best_sol_val = cur_sol_val
                        best_sol = copy(cur_sol)
                    end
                end
            end
            T = cooling_fac * T
            if (time() - start_time) > time_limit  - (n/3)
                break
            end
        end
        println("Current time used: ", time() - start_time)
        println("Currently best solution value: ", best_sol_val)
        
        if best_sol_val < best_sol_bef_MIP_val
            best_sol_bef_MIP_val = best_sol_val
        end
        if best_sol_val < best_ovr_sol_val * 1.00001
            a_vec = compute_a(best_sol, costs, n)
            b_vec = compute_b(best_sol, n, a_vec, len)
            e_vec = calc_e_values(best_sol, costs, n)
            cut_vec = calc_cut_values(best_sol, e_vec, n)
            c_minus = compute_c_minus(cut_vec, n)
            c_plus = compute_c_plus(cut_vec, n)
            dist = compute_distances(best_sol, len, n)
            println("Entering MIP part")
            best_sol, best_sol_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec = local_search_MIP_part(n,len, costs, best_sol, best_sol_val, false, true, Vector{Float64}(), true, [13], a_vec, b_vec, c_minus, c_plus, dist, cut_vec, accleration, time_limit, start_time)
        end
        println("Current time used: ", time() - start_time)
        if best_sol_val < best_ovr_sol_val
            best_ovr_sol_val = best_sol_val
            best_ovr_sol = copy(best_sol)
            n_wo_impr = 0#useless
        else
            n_wo_impr += 1#useless
        end
        cur_sol = randperm(n)
        println("new_random_starting_solution")
        cur_sol_val = calculate_objective_value_new(cur_sol, len, n, costs)
    end
    a_vec = compute_a(best_ovr_sol, costs, n)
    b_vec = compute_b(best_ovr_sol, n, a_vec, len)
    e_vec = calc_e_values(best_ovr_sol, costs, n)
    cut_vec = calc_cut_values(best_ovr_sol, e_vec, n)
    c_minus = compute_c_minus(cut_vec, n)
    c_plus = compute_c_plus(cut_vec, n)
    dist = compute_distances(best_ovr_sol, len, n)
    best_ovr_sol, best_ovr_sol_val, a_vec, b_vec, c_minus, c_plus, dist, cut_vec = local_search_MIP_part(n,len, costs, best_ovr_sol, best_ovr_sol_val, true, true, Vector{Float64}(), true, [17,19], a_vec, b_vec, c_minus, c_plus, dist, cut_vec, accleration, time_limit, start_time)
    println("Total time used: ", time() - start_time)
    println("Number of SA starts: ", n_SA_starts)
    println("real_objecitve_value: ", calculate_objective_value_new(best_ovr_sol, len, n, costs))
    return best_ovr_sol, best_ovr_sol_val, best_sol_bef_MIP_val-best_ovr_sol_val
end


# ---------------------------------------------------------------------------
# Command-line entry point
# Usage:
#   julia --project=. src/window_matheuristic.jl <instance> <cooling_factor> <p_swap> <inner_runs_factor> <acceleration 0|1> [time_limit] [seed]
# Paper settings: 0.976 0.323 87 1
# ---------------------------------------------------------------------------
cooling_fac = parse(Float64, ARGS[2])

T_min = 0.0001
gain_com_prob = 0.4
prob_interchange = parse(Float64, ARGS[3])
n, len, costs = read_instance_data_dummy_palu(ARGS[1])
n_inner_loop = parse(Int64, ARGS[4]) *n

accleration = parse(Float64, ARGS[5])
accleration = Bool(accleration)

# optional 6th argument: time limit in seconds; 0 (or omitted) uses the time limits from the paper
if length(ARGS) >= 6 && parse(Int64, ARGS[6]) > 0
    time_limit = parse(Int64, ARGS[6])
elseif n <= 300
    time_limit = 3600
elseif n <= 500
    time_limit = 7200
else
    time_limit = 10800
end
if time_limit <= n / 2
    @warn "Time limit ($time_limit s) is at most n/2 = $(n / 2) s: the simulated annealing phase will be skipped."
end
# optional 7th argument: random seed
if length(ARGS) >= 7
    Random.seed!(parse(Int64, ARGS[7]))
end

start_time = time()
sol, sol_val, diff_bef_aft_MIP = MSA(cooling_fac, n_inner_loop, T_min, gain_com_prob, prob_interchange, n, costs, len, accleration, time_limit)
duration_time = time()-start_time
filename = basename(ARGS[1])
println(duration_time)
sol_val_wo_acc = sol_val

println(filename, ", ", cooling_fac, ", ", prob_interchange, ", ", ARGS[4], ", ", sol_val, ", ", diff_bef_aft_MIP, ", ", accleration, ", ", time_limit)
