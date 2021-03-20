using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.CombinatorialMaps
using Catlab, Catlab.CategoricalAlgebra.CSets
using CombinatorialSpaces.CombinatorialMaps: TheoryCombinatorialMap
using Catlab.Permutations: cycles

export AbstractPlanarGraph, PlanarGraph, 
    face, src, dst, boundary, isboundary, 
    rem_half_edge!, rem_half_edges!, rem_edge!, 
    add_corolla!,
    planar_rz
    
export incident, σ, α, ϕ, trace_faces, trace_vertices

@present TheoryPlanarGraph <: TheoryCombinatorialMap begin
    V::Ob
    F::Ob

    src::Hom(H,V)   # the source vertex of a half edge
    face::Hom(H,F)  # the face corresponding to a half edge
end

const AbstractPlanarGraph = AbstractACSetType(TheoryPlanarGraph)
const PlanarGraph = CSetType(TheoryPlanarGraph, index=[:ϕ, :src, :face])

face(g::AbstractPlanarGraph, he) = subpart(g, he, :face)
src(g::AbstractPlanarGraph, he) = subpart(g, he, :src)
dst(g::AbstractPlanarGraph, he) = src(g, α(g, he))

"""
    boundary(g::AbstractPlanarGraph)

Returns all half edges on the boundary of a planar graph `g`.
"""
function boundary(g::AbstractPlanarGraph) 
    for f in trace_faces(g)
        if face(g, f[1]) == 0
            return f
        end
    end
    error("Boundary not found")
end

isboundary(g::AbstractPlanarGraph, he) = (face(g, he) == 0)

rem_half_edge!(g::AbstractPlanarGraph, he) = rem_part!(g, :H, he)
rem_half_edges!(g::AbstractPlanarGraph, hes) = rem_parts!(g, :H, hes)

function rem_edge!(g::AbstractPlanarGraph, he)
    
end

function add_corolla!(g::AbstractPlanarGraph, valence::Int)
    v = add_vertex!(g)
    n = nparts(g, :H)
    add_parts!(g, :H, valence; src=v, σ=circshift((n+1):(n+valence), -1))
    return g
end

function pair_half_edges!(g::AbstractPlanarGraph, h1, h2)
    set_subpart!(g, [h1;h2], :α, [h2; h1])
    return g
end

function update_phi!(g::AbstractPlanarGraph)
    hes = collect(1:nparts(g, :H))
    for he in hes
        prev_he = α(g, σ(g, he))
        set_subpart!(g, prev_he, :ϕ, he)
    end
    return g
end

function update_face!(g::AbstractPlanarGraph, he; isboundary = false)
    current_he = he
    next_he = ϕ(g, current_he)
    f = 0
    if !isboundary
        if face(g, he) == 0
            f = add_part!(g, :F)
        else
            f = face(g, he)
        end
    end
    while next_he != he
        set_subpart!(g, current_he, :face, f)
        current_he = next_he
        next_he = ϕ(g, current_he)
    end
    set_subpart!(g, current_he, :face, f)
    return g
end

function planar_rz()
    g = PlanarGraph()
    valences = [3,3,3,3,3,3,3,3,2,4,2,2,2]
    for vl in valences
        add_corolla!(g, vl)
    end
    he_pairs = [(1,26), (2,27), (3,4), (5, 30), (6,7), (8,31), (9,10), (11,33), (12,35),
        (13,36), (14,34), (15,16), (17,32), (18,19), (20,29), (21,22), (23,28), (24,25)]
    pair_half_edges!(g, [p[1] for p in he_pairs], [p[2] for p in he_pairs])
    update_phi!(g)
    update_face!(g, 3; isboundary = true)
    he2f = [2, 5, 8, 11, 20, 1]
    for he in he2f
        update_face!(g, he)
    end

    return g
end
