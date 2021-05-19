using InteractiveUtils


"""
    build_prefix_and_divergence_arrays(H::Array{Int8, 2})

Compute the prefix and divergence arrays for the Positional Burrows-Wheeler Transform on the haplotype matrix H (n_samples x n_snps).
"""
function build_prefix_and_divergence_arrays(H::Array{Int8, 2})
    M, N = size(H)
    u = v = p = q = index = match_start = one(Int32)
    ppa_buffer = Vector{Int32}(undef, M)
    div_buffer = Vector{Int32}(undef, M)
    ppa = Matrix{Int32}(undef, M, N)
    div = Matrix{Int32}(undef, M, N)
    ppa[:, 1] = 1:1:M
    div[:, 1] = ones(Int32, M)
    @inbounds for j in 1:N-1
        u = v = 1
        p = q = j+1
        for i in 1:M
            index = ppa[i, j]
            match_start = div[i, j]
            if match_start > p
                p = match_start
            end
            if match_start > q
                q = match_start
            end
            if H[index, j] == 0
                ppa[u, j+1] = index
                div[u, j+1] = p
                u += 1
                p = 0
            else
                ppa_buffer[v] = index
                div_buffer[v] = q
                v += 1
                q = 0
            end
        end
        for i in u:M
            ppa[i, j+1] = ppa_buffer[i-u+1]
            div[i, j+1] = div_buffer[i-u+1]
        end
    end
    return ppa, div
end
