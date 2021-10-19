

function next(M,i)
    if sum(M[i,:]) == 1
        j = findall(x->x==1,M[i,:])
    else
        j = -1
    end
    return j[1]
end


function cycles(s, e, M)

    N = size(M)[1]
    is_visited = falses(N)
    is_visited[s] = true
    L = []
    l = [s]
    j = next(M,s)
    while j != -1 && !(j in l)
        append!(l,j)
        is_visited[j] = true
        j = next(M,j)
    end

    for i in 1:N
        if is_visited[i] == false && sum(M[i,:]) == 1
            l = [i]
            j = next(M,i)
            while j != -1 && !(j in l)
                append!(l,j)
                is_visited[j] = true
                j = next(M,j)
            end
            push!(L,l)
        end
    end

    return L
end

