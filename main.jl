using JuMP, Gurobi, Distances, LinearAlgebra, Dates
include("read.jl")
include("fonctions.jl")

# Get the instance
inst = readInstance("instance_70_1.txt")
n = inst[1]
d = inst[2]
f = inst[3]
Amin = inst[4]
Nr = inst[5]
R = inst[6]
regions = inst[7]
coords = inst[8]

#Compute distance matrix dij
dij = pairwise(Euclidean(), transpose(coords), transpose(coords))
dij = floor.(Int, dij)

# CHOOSE MODEL TYPE ("poly" or "exp_cycle" or "exp_1" or "exp_2")
model_type = "poly"

####################
# Model definition #
####################

#Model declaration
model = Model(Gurobi.Optimizer)

#Variables declaration
@variable(model, x[i=1:n,j=1:n], Bin) # 1 if the edge (i-j) is in the path

#Objective function
@objective(model, Min, sum(dij .* x)) # Total length of the path

##########################    Constraints   #######################################################
###################################################################################################

@constraint(model, [i in 1:n], x[i,i] == 0)  
@constraint(model, sum(x) >= Amin-1)  #At least Amin airports on the path
@constraint(model, x[dij .> R] .== 0) #  R is the maximal range a plane can fly without landing
@constraint(model, [k in 1:Nr], sum(sum(x[i,j] + x[j,i] for j in 1:n) for i in regions[k]) >= 1) # The path includes airports of every region

@constraint(model, sum(x[d,:]) == 1)
@constraint(model, sum(x[:,d]) == 0)
@constraint(model, sum(x[f,:]) == 0)
@constraint(model, sum(x[:,f]) == 1)

# For each airport on the path there are : 
#    - 1 edge "entering" the airport (except for departure)
#    - 1 edge "leaving" the airport (except for arrival)
@constraint(model, sum(x[d,:]) - sum(x[:,d]) == 1)
@constraint(model, sum(x[f,:]) - sum(x[:,f]) == -1)
@constraint(model, [i in filter(x->!(x in [d,f]), 1:n)], sum(x[i,:]) - sum(x[:,i]) == 0)

###############################################################
#Resolution : 
###############################################################

debut = Dates.Time(Dates.now())

# Nombre de contraintes polynomial
if (model_type == "poly")
    @variable(model, u[1:n], Int)
    @constraint(model, [i in filter(x->!(x in [d]), 1:n), j in filter(x->!(x in [d]), 1:n)], u[j] >= u[i] + 1 + n*(x[i,j]-1))
    JuMP.optimize!(model)
end

# Nombre de contraintes exponentiel
# Sous-problème : recherche de cycles dans un graphe
if (model_type == "exp_cycle")
    while true
        oldstd = stdout
        redirect_stdout(open("null","w"))
        JuMP.optimize!(model)
        redirect_stdout(oldstd)
        println("\n\nCurrent objective value : ", JuMP.objective_value(model))
        cycle=cycles(d,f,value.(x))
        if cycle == []  #Si la fonction cycles ne renvoie rien  i.e pas de sous tour detecté
            return value.(x)
        else 
            for S in cycle
                @constraint(model,sum(sum(x[i,j] for j in S) for i in S) <= size(S,1)-1)
            end
        end
    end
end

# Nombre de contraintes exponentiel
# Sous-problème : version classique des inégalités de sous-tours
if (model_type == "exp_1")
    @variable(model, y[1:n], Bin)
    @constraint(model, [i in filter(x->!(x in [f]), 1:n)], sum(x[i,:]) == y[i])
    @constraint(model, y[f] == 1)

    while true
        oldstd = stdout
        redirect_stdout(open("null","w"))
        JuMP.optimize!(model) # problème maître
        redirect_stdout(oldstd)
        println("\n\nCurrent objective value : ", JuMP.objective_value(model))

        #Model declaration
        model_sep = Model(Gurobi.Optimizer)
        #Variables
        @variable(model_sep, a[i=1:n], Bin)
        #Objective function
        @objective(model_sep, Max, dot(a,value.(x),a) - sum(a) + 1)
        #Constraints
        @constraint(model_sep, sum(a) >= 1)
        #Optimization
        oldstd = stdout
        redirect_stdout(open("null","w"))
        JuMP.optimize!(model_sep)
        redirect_stdout(oldstd)
        
        #Add constraint in master problem
        S = BitArray(round.(value.(a)))
        cardS = sum(round.(value.(a)))
        @constraint(model, sum(x[S,S]) <= cardS - 1)

        if JuMP.objective_value(model_sep) <= 0 # stop condition
            println("\n*********** END OF OPTIMIZATION *************\n")
            return 0
        end
    end
end

# Nombre de contraintes exponentiel
# Sous-problème : version améliorée des inégalités de sous-tours
if (model_type == "exp_2")
    @variable(model, y[1:n], Bin)
    @constraint(model, [i in filter(x->!(x in [f]), 1:n)], sum(x[i,:]) == y[i])
    @constraint(model, y[f] == sum(x[:,f]))

    for v in all_variables(model) # relaxation
        unset_binary(v)
        lb = has_lower_bound(v) ? lower_bound(v) : -Inf
        ub = has_upper_bound(v) ? upper_bound(v) : Inf
        set_lower_bound(v, max(0.0, lb))
        set_upper_bound(v, min(1.0, ub))
    end

    bin = false

    while true
        oldstd = stdout
        redirect_stdout(open("null","w"))
        JuMP.optimize!(model) # problème maître
        redirect_stdout(oldstd)
        println("\n\nCurrent objective value : ", JuMP.objective_value(model))

        #Model declaration
        model_sep = Model(Gurobi.Optimizer)
        #Variables
        @variable(model_sep, a[i=1:n], Bin)
        @variable(model_sep, h[i=1:n], Bin)
        #Objective function
        @objective(model_sep, Max, dot(a,value.(x),a) - dot(value.(y),a) + dot(value.(y),h))
        #Constraints
        @constraint(model_sep, h .<= a)
        @constraint(model_sep, sum(h) == 1)
        #Optimization
        oldstd = stdout
        redirect_stdout(open("null","w"))
        JuMP.optimize!(model_sep)
        redirect_stdout(oldstd)

        #Add constraint in master problem
        S = BitArray(round.(value.(a)))
        i0 = BitArray(round.(value.(h)))
        @constraint(model, sum(x[S,S]) <= sum(y[S]) - sum(y[i0]))

        if bin == true && JuMP.objective_value(model_sep) <= 0 # stop condition
            println("\n*********** END OF OPTIMIZATION *************\n")
            return 0 
        end

        if bin == false && JuMP.objective_value(model_sep) <= 0 # stop relaxation, go back to integer programming
            for v in all_variables(model)
                set_binary(v)
            end 
            global bin = true
        end
    end
end

#Print results
fin = Dates.Time(Dates.now())
obj_value = JuMP.objective_value(model)
println("\n\nObjective value : ", obj_value)
println("Elapsed time : ", Dates.canonicalize(fin - debut))
println("\nVariable x :" )
value.(x)


