using VCFTools


"""
    build_prefix_and_divergence_arrays(H::Array{Int8, 2})

Compute the prefix and divergence arrays (n_samples x n_snps+1) for the 
Positional Burrows-Wheeler Transform on the haplotype matrix H (n_samples x n_snps).
"""
function build_prefix_and_divergence_arrays(H::Array{Int8, 2})
    n_samples, n_snps = size(H)
    u = v = p = q = index = match_start = one(Int32)
    ppa_buffer = Vector{Int32}(undef, n_samples)
    div_buffer = Vector{Int32}(undef, n_samples)
    ppa = Matrix{Int32}(undef, n_samples, n_snps + 1)
    div = Matrix{Int32}(undef, n_samples, n_snps + 1)
    ppa[:, 1] = 1:1:n_samples
    div[:, 1] = ones(Int32, n_samples)
    @inbounds for snp_index in 1:n_snps
        u = v = 1
        p = q = snp_index + 1
        for sample_index in 1:n_samples
            index = ppa[sample_index, snp_index]
            match_start = div[sample_index, snp_index]
            if match_start > p
                p = match_start
            end
            if match_start > q
                q = match_start
            end
            if H[index, snp_index] == 0
                ppa[u, snp_index + 1] = index
                div[u, snp_index + 1] = p
                u += 1
                p = 0
            else
                ppa_buffer[v] = index
                div_buffer[v] = q
                v += 1
                q = 0
            end
        end
        for sample_index in u:n_samples
            ppa[sample_index, snp_index + 1] = ppa_buffer[sample_index - u + 1]
            div[sample_index, snp_index + 1] = div_buffer[sample_index - u + 1]
        end
    end
    return ppa, div
end

struct Segment{T}
    snp_start::T
    snp_end::T
    sample_start::T
    sample_end::T
end

function build_reverse_prefix_array(ppa::Array{Int32, 2})
    M, N = size(ppa)
    reverse_ppa = similar(ppa)
    Threads.@threads for j in 1:N
        @inbounds for i in 1:M
            reverse_ppa[ppa[i, j], j] = i
        end
    end
    return reverse_ppa
end

function build_segments(div::Array{Int32, 2}, reverse_ppa::Array{Int32, 2})
    n_samples = size(div, 1)
    n_snps = size(div, 2) - 1
    segments_by_thread = [Vector{Segment{Int32}}() for t in 1:min(Threads.nthreads(), n_snps)]
    Threads.@threads for snp_index in 1:n_snps
        segments = segments_by_thread[Threads.threadid()]
        last_segment = Segment{Int32}(snp_index, snp_index, convert(Int32, n_snps + 1), convert(Int32, n_snps + 1))
        split_index = convert(Int32, n_snps + 1)
        @inbounds for sample_index in n_samples:-1:1
            match_start = div[sample_index, snp_index + 1]
            if last_segment.snp_start == match_start && last_segment.sample_start < sample_index < last_segment.sample_end
                continue
            elseif match_start == snp_index + 1
                if length(segments) == 0
                    segment = Segment{Int32}(snp_index, snp_index, convert(Int32, sample_index), convert(Int32, sample_index))
                    push!(segments, segment)
                    last_segment = segment
                else
                    for segment_index in length(segments):-1:1
                        segment = segments[segment_index]
                        if segment.sample_start < split_index && segment.snp_start == snp_index
                            break
                        elseif segment.sample_start >= split_index || segment_index == 1
                            segment = Segment{Int32}(snp_index, snp_index, convert(Int32, sample_index), convert(Int32, sample_index))
                            push!(segments, segment)
                            last_segment = segment
                            break
                        end
                    end
                end
                split_index = sample_index
            elseif snp_index == n_snps || match_start != div[reverse_ppa[ppa[sample_index, snp_index + 1], snp_index + 2], snp_index + 2]
                for sample_start in sample_index-1:-1:1
                    if match_start < div[sample_start, snp_index + 1]
                        segment = Segment{Int32}(match_start, snp_index, convert(Int32, sample_start), convert(Int32, sample_index))
                        push!(segments, segment)
                        last_segment = segment
                        break
                    end
                end
            end
        end
    end
    return cat(segments_by_thread..., dims=1)
end

using BenchmarkTools

# H =  Array{Int8, 2}([0 1 0 1 0 1; 1 1 0 0 0 1; 1 1 1 1 1 1; 0 1 1 1 1 0; 0 0 0 0 0 0; 1 0 0 0 1 0; 1 1 0 0 0 1; 0 1 0 1 1 0])
# H =  Array{Int8, 2}([1 0 0 0; 0 0 1 0; 0 0 1 0; 1 0 1 0])
path = "/media/storage/1000_genomes/GRCh38/variants/chr20/yri.chr20.GRCh38.vcf"
H = convert_ht(Int8, path)

ppa, div = build_prefix_and_divergence_arrays(H)
reverse_ppa = build_reverse_prefix_array(ppa)
segments = build_segments(div, reverse_ppa)


H
for i in 1:length(segments)
    segment = segments[i]
    println(segment, " ", H[ppa[segment.sample_end, segment.snp_end + 1], segment.snp_start:segment.snp_end])
end
