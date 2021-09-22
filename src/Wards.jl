data = [
    8.9 14.0 4.3 19.9 2.1 28.0 3.6 1.3 4.3;
    13.5 9.3 4.1 17.5 4.5 26.6 5.7 2.1 4.0;
    18.0 9.9 3.3 19.5 5.7 28.1 4.8 2.4 6.5;
    13.9 10.0 4.7 25.8 2.2 24.0 6.2 1.6 2.9;
    9.5 13.6 3.6 23.4 2.5 22.4 4.2 1.8 3.7;
    13.1 10.1 3.1 23.8 2.3 25.6 2.8 2.4 4.9;
    17.4 5.7 4.7 20.6 4.3 24.3 4.7 3.4 3.3;
    11.4 12.5 4.1 18.8 3.4 18.6 5.2 1.5 3.8
]

N = size(data)[1]
prototypes = copy(data)
clusters = [[i] for i = 1:N]
C = length(clusters)

d²(x⃗, y⃗) = sum((x⃗ - y⃗).^2)
d′(x⃗ᵢ, x⃗ⱼ, nᵢ, nⱼ) = ((nᵢ * nⱼ) / (nᵢ + nⱼ)) * sum((x⃗ᵢ - x⃗ⱼ).^2)

T = 2 # C

for t = 1:T
    D = [i<j ? d′(prototypes[i, :], prototypes[j, :], length(clusters[i]), length(clusters[j])) : (i == j ? Inf : Inf) for i = 1:C, j = 1:C]
    _, idx = findmin(D)
    i, j = idx[1], idx[2]

    append!(clusters[i], clusters[j])
    deleteat!(clusters, j)
    C = length(clusters)

    prototypes = Matrix{Float64}(undef, C, 9)
    for c = 1:C
        vectors = [data[k, :] for k in clusters[c]]
        mean = sum(vectors) / length(vectors)
        prototypes[c, :] = mean
    end

    variances = zeros(C)
    for r = 1:C
        for k in clusters[r]
            variances[r] += d²(data[k, :], prototypes[r, :])
        end
    end
    total_variance = sum(variances)
    println("total variance = $total_variance")
end
