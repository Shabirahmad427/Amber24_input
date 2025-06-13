using PDBTools, MolSimToolkit, Printf
using ProgressMeter: Progress, next!
using StaticArrays: SVector
using LinearAlgebra: norm, dot, diag
using Base.Threads: @threads

export load_simulation, protein_hbonding

function load_simulation(pdbname::String, trjname::String; first=1, step=1, last=nothing)
    sim = MolSimToolkit.Simulation(pdbname, trjname, first=first, step=step, last=last)
    return sim
end

function protein_hbonding(
    simulation::MolSimToolkit.Simulation; ligand="BGLC", hasPBC=false, rOO=3.5, rHO=2.5, α=30.0
)
    monitored = PDBTools.select(
        simulation.atoms, at -> at.resname == ligand && PDBTools.element(at) == "O"
    )
    reference = PDBTools.select(
        PDBTools.select(simulation.atoms, "protein"), at -> !in(PDBTools.element(at), ["H", "C"])
    )
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)
    monitored_hydrogens = filtering_sugar_hydrogens(simulation.atoms, ligand)
    reference_hydrogens = filtering_protein_hydrogens(simulation.atoms)
    XYZbuffer = Matrix{Float32}(undef, length(simulation.atoms), 3)
    hbonds = Dict{Tuple{String, String}, Int64}()
    uc = !hasPBC ? guess_unitcell(simulation.atoms) : nothing
    nframes = length(simulation.frame_range)
    progress = Progress(nframes)
    for frame in simulation
        xyz = MolSimToolkit.positions(frame)
        uc = hasPBC ? diag(MolSimToolkit.unitcell(frame)) : uc
        # 1st H-bond criteria trial
        hbond_1st_criteria = MolSimToolkit.minimum_distances(
            xpositions = [ SVector(xyz[i]) for i in imonitored ],
            ypositions = [ SVector(xyz[j]) for j in jreference ],
            xn_atoms_per_molecule = 1,
            cutoff = rOO,
            unitcell = uc
        )
        candidates = findall(md -> md.within_cutoff, hbond_1st_criteria)
        isempty(candidates) && continue
        @inbounds for i in eachindex(xyz)
            XYZbuffer[i, 1] = xyz[i][1]
            XYZbuffer[i, 2] = xyz[i][2]
            XYZbuffer[i, 3] = xyz[i][3]
        end
        # moving the candidates to the final checking
        hbonds_lock = ReentrantLock()
        @inbounds @threads :static for k in eachindex(candidates)
            ligand_atoms = candidates[k]
            i = imonitored[hbond_1st_criteria[ligand_atoms].i]  # the expected substrate
            j = jreference[hbond_1st_criteria[ligand_atoms].j]  # and protein atoms bonded by this interaction
            oa = @view(XYZbuffer[j,:])
            od = MolSimToolkit.wrap(@view(XYZbuffer[i,:]), oa, uc)
            if haskey(monitored_hydrogens, i)   ## ligand as h-bond donor
                hydrogen = monitored_hydrogens[i]
                hd = MolSimToolkit.wrap(@view(XYZbuffer[hydrogen,:]), oa, uc)
                if hbond_2nd_criteria(hd, oa, rHO) && hbond_3rd_criteria(hd, od, oa, α)
                    hdonor_residue = join([PDBTools.resname(simulation.atoms[hydrogen]), PDBTools.resnum(simulation.atoms[hydrogen])], "")
                    labeld = join([hdonor_residue, PDBTools.name(simulation.atoms[hydrogen])], "-")
                    haccep_residue = join([PDBTools.resname(simulation.atoms[j]), PDBTools.resnum(simulation.atoms[j])], "")
                    labela = join([haccep_residue, PDBTools.name(simulation.atoms[j])], "-")
                    lock(hbonds_lock) do
                        if haskey(hbonds, (labeld, labela))
                            hbonds[(labeld, labela)] += 1
                        else
                            hbonds[(labeld, labela)] = 1
                        end
                    end
                end
            end
            if haskey(reference_hydrogens, j)   ## ligand as h-bond acceptor
                hydrogens = reference_hydrogens[j]
                @fastmath for hydrogen in hydrogens
                    hd = MolSimToolkit.wrap(@view(XYZbuffer[hydrogen, :]), od, uc)
                    if hbond_2nd_criteria(hd, od, rHO) && hbond_3rd_criteria(hd, oa, od, α)
                        hdonor_residue = join([PDBTools.resname(simulation.atoms[hydrogen]), PDBTools.resnum(simulation.atoms[hydrogen])], "")
                        labeld = join([hdonor_residue, PDBTools.name(simulation.atoms[hydrogen])], "-")
                        haccep_residue = join([PDBTools.resname(simulation.atoms[i]), PDBTools.resnum(simulation.atoms[i])], "")
                        labela = join([haccep_residue, PDBTools.name(simulation.atoms[i])], "-")
                        lock(hbonds_lock) do
                            if haskey(hbonds, (labeld, labela))
                                hbonds[(labeld, labela)] += 1
                            else
                                hbonds[(labeld, labela)] = 1
                            end
                        end
                    end
                end
            end
        end
        next!(progress)
    end

    Base.open("./output_hbonding.dat", "w") do file
        println(file, "Found $(length(hbonds)) hbonds.")
        println(file, "Donor       Acceptor   |  Occupancy")
        for (key, value) in hbonds
            info = Printf.@sprintf("%15s %15s  |    %5.2f", key[1], key[2], 100 * value / nframes)
            println(info)
            println(file, info)
        end
    end
    return nothing
end

@inline function filtering_sugar_hydrogens(atoms::Vector{<:PDBTools.Atom}, ligand::String; tol::Float64=1.0)
    bonds = Dict{Int32, Int32}()
    O = PDBTools.select(atoms, at -> at.resname == ligand && PDBTools.element(at) == "O")
    H = PDBTools.select(atoms, at -> at.resname == ligand && PDBTools.element(at) == "H")
    for Oi in O
        Hi = PDBTools.select(H, at -> at.resname == Oi.resname && norm(PDBTools.coor(at) - PDBTools.coor(Oi)) <= tol)
        idx = PDBTools.index.(Hi)
        if !isempty(idx)
            bonds[Oi.index] = Int32(idx[1])
        end
    end
    return bonds
end

@inline function filtering_protein_hydrogens(atoms::Vector{<:PDBTools.Atom}; tol::Float64=1.5)
    bonds = Dict{Int32, Vector{Int32}}()
    protein_atoms = PDBTools.select(atoms, "protein")
    X = PDBTools.select(protein_atoms, at -> !in(PDBTools.element(at), ["H", "C"]))
    H = PDBTools.select(protein_atoms, at -> PDBTools.element(at) == "H")
    for Xi in X
        Hi = PDBTools.select(H, at -> at.resname == Xi.resname && norm(PDBTools.coor(at) - PDBTools.coor(Xi)) <= tol)
        idx = PDBTools.index.(Hi)
        if !isempty(idx)
            bonds[Xi.index] = Int32.(idx)
        end
    end
    return bonds
end

@inline function guess_unitcell(atoms::Vector{<:PDBTools.Atom}; padding::Float64=0.0)
    return PDBTools.maxmin(atoms).xlength .+ padding
end

@inline function hbond_2nd_criteria(HD::AbstractVector, OA::AbstractVector, cutoff::Float64)
    v1 = HD .- OA
    return norm(v1) <= cutoff
end

@inline function hbond_3rd_criteria(HD::AbstractVector, OD::AbstractVector, OA::AbstractVector, cutoff::Float64)
    v1 = HD .- OD
    v2 = OA .- OD
    α = hbond_angle(v1, v2)
    return α <= cutoff
end

@inline function hbond_angle(v1::AbstractVector, v2::AbstractVector)
    ratio = dot(v1, v2) / (norm(v1) * norm(v2))
    ratio = clamp(ratio, -1.0, 1.0)
    return π \ 180 * acos(ratio)
end

# function ligand_hbonding(
#     simulation::MolSimToolkit.Simulation; selection="not water", water="TIP3", rOO=3.5, rHO=2.5, α=30.0
# )
#     monitored = PDBTools.select(simulation.atoms, at -> at.resname == water && PDBTools.element(at) != "H")
#     reference = PDBTools.select(
#         PDBTools.select(simulation.atoms, selection),
#         at -> PDBTools.element(at) != "H"
#     )
#     imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)
#     monitored_hydrogens = filtering_ligand_hydrogens(simulation.atoms, imonitored, water)
#     reference_hydrogens = filtering_sugar_hydrogens(simulation.atoms, jreference, water)
#     M = falses(length(imonitored), length(simulation.frame_range))
#     XYZbuffer = Matrix{Float32}(undef, length(simulation.atoms), 3)
#     for (iframe, frame) in enumerate(simulation)
#         xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
#         mindist = MolSimToolkit.minimum_distances(
#             xpositions = [ SVector(xyz[i]) for i in imonitored ],
#             ypositions = [ SVector(xyz[j]) for j in jreference ],
#             xn_atoms_per_molecule = 1,
#             cutoff = rOO,
#             unitcell = uc
#         )
#         candidates = findall(md -> md.within_cutoff, mindist)
#         isempty(candidates) && continue
#         @inbounds for i in eachindex(xyz)
#             XYZbuffer[i, 1] = xyz[i][1]
#             XYZbuffer[i, 2] = xyz[i][2]
#             XYZbuffer[i, 3] = xyz[i][3]
#         end
#         @inbounds @threads :static for k in eachindex(candidates)
#             iwater = candidates[k]
#             i, j = imonitored[mindist[iwater].i], jreference[mindist[iwater].j]
#             oa = @view(XYZbuffer[j, :])
#             od = MolSimToolkit.wrap(@view(XYZbuffer[i, :]), oa, uc)
#             if haskey(monitored_hydrogens, i)   ## water as h-bond donor
#                 hbond_found = false
#                 @fastmath for hydrogen in monitored_hydrogens[i]
#                     hd = MolSimToolkit.wrap(@view(XYZbuffer[hydrogen, :]), oa, uc)
#                     hbond_found = hbond_bond2_checking(hd, oa, rHO) && hbond_angle_checking(hd, od, oa, α)
#                     hbond_found && break
#                 end
#                 M[iwater, iframe] = hbond_found
#                 hbond_found && continue
#             end
#             if haskey(reference_hydrogens, j)   ## water as h-bond acceptor
#                 hydrogen = reference_hydrogens[j]
#                 hd = MolSimToolkit.wrap(@view(XYZbuffer[hydrogen, :]), od, uc)
#                 M[iwater, iframe] = hbond_bond2_checking(hd, od, rHO) && hbond_angle_checking(hd, oa, od, α)
#             end
#         end
#     end
#     return BitMatrix(M)
# end
