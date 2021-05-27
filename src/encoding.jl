using DataStructures
using VCFTools


struct Segment{T}
    snp_start::T
    snp_end::T
    sample_start::T
    sample_end::T
end

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

function build_segments(div::Array{Int32, 2}, ppa::Array{Int32, 2}, reverse_ppa::Array{Int32, 2})
    n_samples = size(div, 1)
    n_snps = size(div, 2) - 1
    segments_by_thread = [Vector{Segment{Int32}}() for t in 1:min(Threads.nthreads(), n_snps)]
    Threads.@threads for snp_index in 1:n_snps
        segments = segments_by_thread[Threads.threadid()]
        last_segment = Segment{Int32}(snp_index, snp_index, convert(Int32, n_samples + 1), convert(Int32, n_samples + 1))
        split_index = convert(Int32, n_samples + 1)
        @inbounds for sample_index in n_samples:-1:1
            match_start = div[sample_index, snp_index + 1]
            if match_start == snp_index
                continue
            elseif last_segment.snp_start == match_start && last_segment.sample_start < sample_index < last_segment.sample_end
                continue
            elseif match_start == snp_index + 1
                segment = Segment{Int32}(snp_index, snp_index, convert(Int32, sample_index), convert(Int32, split_index - 1))
                push!(segments, segment)
                last_segment = segment
                split_index = sample_index
            elseif snp_index == n_snps || match_start != div[reverse_ppa[ppa[sample_index, snp_index + 1], snp_index + 2], snp_index + 2]
                sample_start = sample_index
                sample_end = n_samples
                for neighbor_index in sample_index-1:-1:1
                    if match_start < div[neighbor_index, snp_index + 1]
                        sample_start = neighbor_index
                        break
                    end
                end
                for neighbor_index in sample_index+1:n_samples
                    if match_start < div[neighbor_index, snp_index + 1]
                        sample_end = neighbor_index - 1
                        break
                    end
                end
                segment = Segment{Int32}(match_start, snp_index, convert(Int32, sample_start), convert(Int32, sample_end))
                push!(segments, segment)
                last_segment = segment
            end
        end
    end
    return cat(segments_by_thread..., dims=1)
end

function build_segment_overlaps(segments::Vector{Segment{Int32}}, ppa::Array{Int32, 2})
    n_segments = length(segments)
    n_snps = size(ppa, 2) - 1
    segment_index_by_snp = [Vector{Int32}() for snp_index in 1:n_snps]
    segment_overlaps = [Vector{Int32}() for segment_index in 1:n_segments]
    for segment_index in 1:n_segments
        segment = segments[segment_index]
        for snp_index in segment.snp_start:segment.snp_end
            push!(segment_index_by_snp[snp_index], segment_index)
        end
    end
    Threads.@threads for segment_index in 1:n_segments
        segment = segments[segment_index]
        candidate_overlap_set = Set{Int32}()
        original_form_sample_indices = Vector{Int32}()
        @inbounds for sample_index in segment.sample_start:segment.sample_end
            push!(original_form_sample_indices, ppa[sample_index, segment.snp_end + 1])
        end
        @inbounds for snp_index in segment.snp_start:segment.snp_end
            for candidate_index in segment_index_by_snp[snp_index]
                # if segment_index == candidate_index || in(candidate_index, candidate_overlap_set)
                if in(candidate_index, candidate_overlap_set)
                    continue
                end
                push!(candidate_overlap_set, candidate_index)
                candidate_segment = segments[candidate_index]
                for sample_index in candidate_segment.sample_start:candidate_segment.sample_end
                    if in(ppa[sample_index, candidate_segment.snp_end + 1], original_form_sample_indices)
                        push!(segment_overlaps[segment_index], candidate_index)
                        break
                    end
                end
            end
        end
    end
    return segment_overlaps
end

function calculate_scores(coverage::Array{Int8, 2}, ppa::Array{Int32, 2}, segments::Vector{Segment{Int32}}, indices::Vector{Int32})
    scores = Vector{Int32}(undef, length(indices))
    Threads.@threads for (score_index, segment_index) in collect(enumerate(indices))
        segment = segments[segment_index]
        @inbounds original_form_sample_indices = [ppa[sample_index, segment.snp_end + 1] for sample_index in segment.sample_start:segment.sample_end]
        score = zero(Int32)
        for snp_index in segment.snp_start:segment.snp_end
            for original_form_sample_index in original_form_sample_indices
                @inbounds score += coverage[original_form_sample_index, snp_index]
            end
        end
        @inbounds scores[score_index] = score
    end
    return scores
end

function calculate_scores(coverage::Array{Int8, 2}, ppa::Array{Int32, 2}, segments::Vector{Segment{Int32}})
    indices = convert(Vector{Int32}, 1:length(segments) |> collect)
    return calculate_scores(coverage, ppa, segments, indices)
end

function update_coverage!(coverage::Array{Int8, 2}, ppa::Array{Int32, 2}, segments::Vector{Segment{Int32}}, segment_index::Int)
    segment = segments[segment_index]
    original_form_sample_indices = [ppa[sample_index, segment.snp_end + 1] for sample_index in segment.sample_start:segment.sample_end]
    for snp_index in segment.snp_start:segment.snp_end
        for original_form_sample_index in original_form_sample_indices
            @inbounds coverage[original_form_sample_index, snp_index] = 0
        end
    end
end

function update_scores!(scores::MutableBinaryMaxHeap{Int32}, coverage::Array{Int8, 2}, ppa::Array{Int32, 2}, segments::Vector{Segment{Int32}}, indices::Vector{Int32})
    updated_scores = calculate_scores(coverage, ppa, segments, indices)
    for (segment_index, score) in zip(indices, updated_scores)
        update!(scores, convert(Int, segment_index), score)
    end
end

function build_encoding(H::Array{Int8, 2}, max_size::Int)
    encoding = Vector{Segment{Int32}}(undef, 0)
    ppa, div = build_prefix_and_divergence_arrays(H)
    reverse_ppa = build_reverse_prefix_array(ppa)
    segments = build_segments(div, ppa, reverse_ppa)
    segment_overlaps = build_segment_overlaps(segments, ppa)
    coverage = ones(Int8, size(H, 1), size(H, 2))
    scores = MutableBinaryMaxHeap{Int32}(calculate_scores(coverage, ppa, segments))
    score, segment_index = top_with_handle(scores)
    while score > 0 && length(encoding) < max_size
        push!(encoding, segments[segment_index])
        update_coverage!(coverage, ppa, segments, segment_index)
        update_indices = [index for index in segment_overlaps[segment_index] if scores[convert(Int, index)] > 0]
        update_scores!(scores, coverage, ppa, segments, update_indices)
        score, segment_index = top_with_handle(scores)
    end
    return encoding, ppa
end

function build_encoding(H::Array{Int8, 2})
    return build_encoding(H, length(H))
end

function build_encoding_haplotypes(H::Array{Int8, 2}, ppa::Array{Int32, 2}, encoding::Vector{Segment{Int32}})
    n_segments = length(encoding)
    encoding_haplotypes = Vector{Vector{Int8}}(undef, n_segments)
    Threads.@threads for segment_index in 1:n_segments
        @inbounds segment = encoding[segment_index]
        @inbounds encoding_haplotypes[segment_index] = H[ppa[segment.sample_start, segment.snp_end + 1], segment.snp_start:segment.snp_end]
    end
    return encoding_haplotypes
end

using BenchmarkTools

H =  Array{Int8, 2}([0 1 0 1 0 1; 1 1 0 0 0 1; 1 1 1 1 1 1; 0 1 1 1 1 0; 0 0 0 0 0 0; 1 0 0 0 1 0; 1 1 0 0 0 1; 0 1 0 1 1 0])
# H =  Array{Int8, 2}([1 0 0 0; 0 0 1 0; 0 0 1 0; 1 0 1 0])
path = "/media/storage/1000_genomes/GRCh38/variants/chr20/yri.chr20.GRCh38.vcf"
H = convert_ht(Int8, path)
H = H[:, 1:100000]
@time encoding, ppa = build_encoding(H)
@time encoding_haplotypes = build_encoding_haplotypes(H, ppa, encoding)


for segment in encoding
    println(segment)
    println(H[ppa[segment.sample_start:segment.sample_end, segment.snp_end + 1], 1:segment.snp_end])
end


ppa, div = build_prefix_and_divergence_arrays(H)
reverse_ppa = build_reverse_prefix_array(ppa)
segments = build_segments(div, ppa, reverse_ppa)
segment_overlaps = build_segment_overlaps(segments, ppa)
coverage = ones(Int8, size(H, 1), size(H, 2))
scores = MutableBinaryMaxHeap{Int32}(calculate_scores(coverage, ppa, segments))

segment_index = 15
segment = segments[segment_index]
H[ppa[:, segment.snp_end + 1], :]
for segment_index in segment_overlaps[segment_index]
    segment = segments[segment_index]
    println(segment)
    println(H[ppa[:, segment.snp_end + 1], 1:segment.snp_end])
end

