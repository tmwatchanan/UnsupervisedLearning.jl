X = Float64[2 5; 6 4; 5 3; 2 2; 1 4; 5 4; 3 3; 2 3; 2 4; 8 2; 9 2; 10 2; 11 2; 10 3; 9 1]
d2(x⃗, y⃗) = sqrt(sum((x⃗ - y⃗).^2))
d2²(x⃗, y⃗) = sum((x⃗ - y⃗).^2)
N = size(X)[1]

# prototypes = X[1:4, :]
prototypes = X[12:15, :]
K = size(prototypes)[1]

function compute_distances(X, prototypes)
    N = size(X)[1]
    K = size(prototypes)[1]

    distances = Matrix{Float64}(undef, N, K)
    for k = 1:K
        for n = 1:N
            distances[n, k] = d2(X[n, :], prototypes[k, :])
        end
    end
    distances
end

for t = 1:10
    previous_prototypes = copy(prototypes)

    distances = compute_distances(X, prototypes)
    
    cluster_assignments = map(argmin, eachrow(distances))
    println("cluster t=$t @ ", cluster_assignments)
    for k = 1:K
        indices = findall(x->x==k, cluster_assignments)
        prototypes[k, :] = sum(X[indices, :], dims=1) / length(indices)
    end
    if all(prototypes == previous_prototypes)
        break
    end
end

distances = compute_distances(X, prototypes)
cluster_assignments = map(argmin, eachrow(distances))
println("cluster @ ", cluster_assignments)

J = 0.0
for k = 1:K
    indices = findall(x->x==k, cluster_assignments)
    for i = 1:length(indices)
        J += d2²(X[indices[i], :], prototypes[k, :])
    end
end
println("J = ", J)
