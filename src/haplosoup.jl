using DataStructures


include("algorithms.jl")

using Distributed
using BenchmarkTools
using Test
using VCFTools

path = "/media/storage/1000_genomes/GRCh38/variants/chr20/yri.chr20.GRCh38.vcf"
H = convert_ht(Int8, path)
H =  Array{Int8, 2}([0 1 0 1 0 1; 1 1 0 0 0 1; 1 1 1 1 1 1; 0 1 1 1 1 0; 0 0 0 0 0 0; 1 0 0 0 1 0; 1 1 0 0 0 1; 0 1 0 1 1 0])


function build_scoring_blocks(div::Array{Int32, 2})
    M, N = size(div)
    block_start = Matrix{Int32}(undef, M, N)
    block_end = Matrix{Int32}(undef, M, N)
    @sync Threads.@threads for j in 1:N
        segment_length_i = segment_length_neighbor = 0
        block_start[1, j] = 1
        block_end[1, j] = 1
        @inbounds for i in 2:M
            segment_length_i = j - div[i, j]
            if segment_length_i == 0
                block_start[i, j] = i
                block_end[i, j] = i
                continue
            end
            neighbor_i = 1
            for outer neighbor_i in i-1:-1:1
                segment_length_neighbor = j - div[neighbor_i, j]
                if segment_length_i >= segment_length_neighbor
                    break
                end
            end
            if segment_length_i == segment_length_neighbor
                block_start[i, j] = block_start[neighbor_i, j]
                block_end[i, j] = block_end[neighbor_i, j]
                continue
            elseif segment_length_i < segment_length_neighbor
                block_start[i, j] = neighbor_i + 1
            else
                block_start[i, j] = neighbor_i
            end
            if i == M
                block_end[i, j] = i
                continue
            end
            for outer neighbor_i in i+1:M
                segment_length_neighbor = j - div[neighbor_i, j]
                if segment_length_i > segment_length_neighbor
                    break
                end
            end
            if segment_length_i > segment_length_neighbor
                block_end[i, j] = neighbor_i - 1
            else
                block_end[i, j] = neighbor_i
            end
        end
    end
    return block_start, block_end
end

function update_scores!(
    scores::MutableBinaryMaxHeap{Int32},
    coverage::Array{Int8, 2},
    div::Array{Int32, 2},
    ppa::Array{Int32, 2},
    reverse_ppa::Array{Int32, 2},
    block_start::Array{Int32, 2},
    block_end::Array{Int32, 2},
    i_range_start::Int,
    i_range_end::Int,
    j_range_start::Int,
    j_range_end::Int,
    j_selected::Int,
    )
    M = size(coverage, 1)
    @inbounds for j in j_range_start:j_range_end
        for i in i_range_start:i_range_end
            score = zero(Int32)
            i_adjusted = reverse_ppa[ppa[i, j_selected + 1], j+1]
            for j_score in div[i_adjusted, j+1]:j
                for i_score in block_start[i_adjusted, j+1]:block_end[i_adjusted, j+1]
                    score += coverage[ppa[i_score, j+1], j_score]
                end
            end
            update!(scores, (j - 1) * M + i_adjusted, score)
        end
    end
    return nothing
end

# function update_scores!(
#     scores::MutableBinaryMaxHeap{Int32},
#     coverage::Array{Int8, 2},
#     div::Array{Int32, 2},
#     ppa::Array{Int32, 2},
#     block_start::Array{Int32, 2},
#     block_end::Array{Int32, 2},
#     i_range_start::Int,
#     i_range_end::Int,
#     j_range_start::Int,
#     j_range_end::Int
#     )
#     for j in j_range_start:j_range_end
#         M = size(coverage, 1)
#         @inbounds for i in i_range_start:i_range_end
#             score = zero(Int32)
#             for block_j in div[i, j+1]:j
#                 for block_i in block_start[i, j+1]:block_end[i, j+1]
#                     score += coverage[ppa[block_i, block_j + 1], block_j]
#                 end
#                 println()
#                 println(i, " ", j)
#                 println(score)
#                 println(sum(coverage[ppa[block_start[i, j+1]:block_end[i, j+1], block_j+1], block_j]))
#             end
#             println("HERE")
#             # println(score)
#             # println(scores[(j - 1) * M + i])
#             update!(scores, (j - 1) * M + i, score)
#             # println(scores[(j - 1) * M + i])
#             # println()
#         end
#     end
#     return nothing
# end

function update_scores!(
    scores::MutableBinaryMaxHeap{Int32},
    coverage::Array{Int8, 2},
    div::Array{Int32, 2},
    ppa::Array{Int32, 2},
    reverse_ppa::Array{Int32, 2},
    block_start::Array{Int32, 2},
    block_end::Array{Int32, 2}
    )
    update_scores!(scores, coverage, div, ppa, reverse_ppa, block_start, block_end, 1, size(coverage, 1), 1, size(coverage, 2), size(coverage, 2))
end

function build_reverse_prefix_array(ppa::Array{Int32, 2})
    M, N = size(ppa)
    reverse_ppa = similar(ppa)
    @sync Threads.@threads for j in 1:N
        @inbounds for i in 1:M
            reverse_ppa[ppa[i, j], j] = i
        end
    end
    return reverse_ppa
end

function find_altered_block(
    div::Array{Int32, 2},
    ppa::Array{Int32, 2},
    reverse_ppa::Array{Int32, 2},
    block_start::Array{Int32, 2},
    block_end::Array{Int32, 2},
    i_selected::Int,
    j_selected::Int
    )
    N = size(div, 2) - 1
    i_start = convert(Int, block_start[i_selected, j_selected + 1])
    i_end = convert(Int, block_end[i_selected, j_selected + 1])
    # j_start = convert(Int, minimum(div[i_start:i_end, j_selected + 1]))
    j_start = convert(Int, div[i_selected, j_selected + 1])
    j_end_candidates = Vector{Int}(undef, i_end - i_start + 1)
    @sync Threads.@threads for i in i_start:i_end
        original_index = ppa[i, j_selected + 1]
        div_i_j_selected = div[i, j_selected + 1]
        j_end_candidates[i - i_start + 1] = j_selected
        j = j_selected
        @inbounds for outer j in j_selected+1:N
            if div[reverse_ppa[original_index, j+1], j+1] > div_i_j_selected
                break
            end
        end
        if div[reverse_ppa[original_index, j+1], j+1] > div_i_j_selected
            j_end_candidates[i - i_start + 1] = j - 1
        else
            j_end_candidates[i - i_start + 1] = j
        end
    end
    j_end = maximum(j_end_candidates)
    return i_start, i_end, j_start, j_end
end

function find_selected_segment(
    H::Array{Int8, 2},
    div::Array{Int32, 2},
    ppa::Array{Int32, 2},
    i_selected::Int,
    j_selected::Int
    )
    return H[ppa[i_selected, j_selected + 1], div[i_selected, j_selected + 1]:j_selected]
end

function update_coverage!(
    coverage::Array{Int8, 2},
    div::Array{Int32, 2},
    ppa::Array{Int32, 2},
    block_start::Array{Int32, 2},
    block_end::Array{Int32, 2},
    i_selected::Int,
    j_selected::Int
    )
    @sync Threads.@threads for j in div[i_selected, j_selected + 1]:j_selected
        @inbounds for i in block_start[i_selected, j_selected + 1]:block_end[i_selected, j_selected + 1]
            coverage[ppa[i, j_selected + 1], j] = 0
        end
    end
end


H = convert_ht(Int8, path)
# H = Hbig[:, 1:10000]
ppa, div = build_prefix_and_divergence_arrays(H)
reverse_ppa = build_reverse_prefix_array(ppa)
block_start, block_end = build_scoring_blocks(div)
coverage = ones(Int8, size(H)...)
scores = MutableBinaryMaxHeap{Int32}(zeros(Int32, length(H)))
@time update_scores!(scores, coverage, div, ppa, reverse_ppa, block_start, block_end)
indices = CartesianIndices(size(H))
score, k = top_with_handle(scores)
# while score > 0
#     i_selected, j_selected = indices[k][1], indices[k][2]
#     selected_segment = find_selected_segment(H, div, ppa, i_selected, j_selected)
#     println(score, " ", div[i_selected, j_selected + 1], " ", selected_segment)
#     i_start, i_end, j_start, j_end = find_altered_block(div, ppa, reverse_ppa, block_start, block_end, i_selected, j_selected)
#     update!(scores, k, 0)
#     update_coverage!(coverage, div, ppa, block_start, block_end, i_selected, j_selected)
#     update_scores!(scores, coverage, div, ppa, reverse_ppa, block_start, block_end, i_start, i_end, j_start, j_end, j_selected)
#     score, k = top_with_handle(scores)
# end

@time begin
    score, k = top_with_handle(scores)
    #
    i_selected, j_selected = indices[k][1], indices[k][2]
    selected_segment = find_selected_segment(H, div, ppa, i_selected, j_selected)
    # println(score, " ", div[i_selected, j_selected + 1], " ", selected_segment)
    i_start, i_end, j_start, j_end = find_altered_block(div, ppa, reverse_ppa, block_start, block_end, i_selected, j_selected)
    update!(scores, k, 0)
    update_coverage!(coverage, div, ppa, block_start, block_end, i_selected, j_selected)
    update_scores!(scores, coverage, div, ppa, reverse_ppa, block_start, block_end, i_start, i_end, j_start, j_end, j_selected)
    # score, k = top_with_handle(scores)
    # coverage
end





H
i_selected, j_selected = indices[k][1], indices[k][2]
selected_segment = find_selected_segment(H, div, ppa, i_selected, j_selected)
i_start, i_end, j_start, j_end = find_altered_block(div, ppa, reverse_ppa, block_start, block_end, i_selected, j_selected)
coverage
update_coverage!(coverage, div, ppa, block_start, block_end, i_selected, j_selected)
coverage
update_scores!(scores, coverage, div, ppa, block_start, block_end, i_start, i_end, j_start, j_end)
scores

top_with_handle(scores)
i, j = indices[35][1], indices[35][2]
H[ppa[:, j+1], 1:j]
div[:, j+1]
coverage[ppa[:, j+1], 1:j]



scores = MutableBinaryMaxHeap{Int32}(zeros(Int32, length(H)))
coverage = ones(Int8, size(H)...)
update_scores!(scores, coverage, div, ppa, block_start, block_end)
i_selected, j_selected = indices[46][1], indices[46][2]
update_coverage!(coverage, div, ppa, block_start, block_end, i_selected, j_selected)
print("START")
H[ppa[:, j_selected + 1], 1:j_selected]
for j in 1:j_selected
    println(coverage[ppa[:, j+1], j])
end
i_start, i_end, j_start, j_end = find_altered_block(div, ppa, reverse_ppa, block_start, block_end, i_selected, j_selected)
j_selected
###
for j in j_start:j_end
    for i in i_start:i_end
        println(ppa[i, j_selected + 1])
        println(reverse_ppa[ppa[i, j_selected + 1], j+1], " ", j)
    end
end
###
reverse_ppa[ppa[6, 7], 6]
coverage[reverse_ppa[ppa[6, 7], 6], 5]
H[ppa[:, 6], 1:5]
#
scores
update_scores!(scores, coverage, div, ppa, block_start, block_end, i_start, i_end, j_start, j_end)
scores
#
for j in j_start:j_end
    # println(H[ppa[i_start:i_end, j+1], j])
    println(coverage[ppa[i_start:i_end, j+1], j])
    # coverage[ppa[i_start:i_end, j+1], j] = 0
    # coverage[ppa[i, j+1], j] = 0
end

score = 0
for j in div[3, 6]:5
    println(coverage[ppa[block_start[3, 6]:block_end[3, 6], j+1], j])
end
coverage[:, 1:5]

H[ppa[6, 7], 1:6]
reverse_ppa[]



H[ppa[:, j_selected+1], 1:j_selected]
find_altered_block(div, ppa, reverse_ppa, block_start, block_end, i_selected, j_selected)

i_selected = convert(Int32, 8)
j_selected = convert(Int32, 2)
H[ppa[:, j_selected+1], 1:j_selected]
find_altered_block(div, ppa, reverse_ppa, block_start, block_end, i_selected, j_selected)
