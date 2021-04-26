module TBhelpers

using LinearAlgebra
using ACEtb.SlaterKoster

export set_SK_orbitals, find_row_in_matrix, get_latvec_ids, get_R, get_tbcells, wrap_shift, get_sparse_indexing, get_neighbours, get_translation_cells

function set_SK_orbitals(ntypes,nshells,orblist)
   orb_keys = Dict(0 => "s", 1 => "p", 2 => "d", 3 => "f")
   # Defining orbital keys for SK package
   k = 0
   skokeys = []
   for i = 1:ntypes
      okeys = []
      for j = 1:nshells[i]
         k += 1
         push!(okeys, SKOrbital(orb_keys[orblist[k]]))
      end
      push!(skokeys, okeys)
   end

   # generate bonds automatically
   SKH_list = []
   for ki in skokeys
       push!(SKH_list, SKH([i for i in ki]))
   end
   return SKH_list
end

function find_row_in_matrix(row, mat; round_digits=nothing)
    if round_digits != nothing
        row = round.(row, digits=round_digits)
    end
    rtn = []
    for i=1:size(mat,1)
        if round_digits != nothing
            matrow = round.(mat[i,:], digits=round_digits)
        else
            matrow = mat[i,:]
        end
        if row == matrow
            push!(rtn,i)
        end
    end
    return rtn
end

function get_latvec_ids(x,y,z)
    # Lattice vectors
    inds = [CartesianIndices((x,y,z))...]
    indx = [convert(Array{Int64,1}, [d[1],d[2],d[3]]) for d in inds]
    R_cN = hcat(indx...)
    N_c = Array{Int64,2}([x y z]')
    R_cN .+= N_c .÷ 2
    R_cN .%= N_c
    R_cN .-= N_c .÷ 2
    return R_cN
end

function get_latvec_ids(v::AbstractArray)
    # Lattice vectors
    inds = [CartesianIndices(tuple(v...))...]
    indx = [convert(Array{Int64,1}, [d[1],d[2],d[3]]) for d in inds]
    R_cN = hcat(indx...)
    N_c = Array{Int64,2}(v')
    R_cN .+= N_c .÷ 2
    R_cN .%= N_c 
    R_cN .-= N_c .÷ 2
    return R_cN'
end

function get_latvec_ids(vmin::AbstractArray, vmax::AbstractArray; origin=true)
    if origin
       R_cN = [[0 0 0]] 
    else
       R_cN = []
    end 
    for i = vmin[1]:vmax[1]
       for j = vmin[1]:vmax[1]
          for k = vmin[1]:vmax[1]
              if ( origin && i==0 && j==0 && k==0 ) 
                 continue 
              end 
              push!(R_cN, [i j k]) 
          end 
       end
    end
    return R_cN
end

function get_latvec_ids(v::AbstractArray)
    @assert length(v) < 3
    # Lattice vectors
    inds = [CartesianIndices(tuple(v))...]
    indx = [convert(Array{Int64,1}, [d[1],d[2],d[3]]) for d in inds]
    R_cN = hcat(indx...)
    N_c = Array{Int64,2}(v')
    R_cN .+= N_c .÷ 2
    R_cN .%= N_c
    R_cN .-= N_c .÷ 2
    return R_cN
end

function get_R(atoms, i, j; cell_vec=nothing)
    # i is always at unitcell (0,0,0)
    # @show atoms
    # @show i,j
    pos = positions(atoms)
    # @show pos
    if cell_vec != nothing
        c = cell(atoms)
        return pos[j] + (c * cell_vec) - pos[i]
    else
        return pos[j] - pos[i]
    end
end

function get_translation_cells(cell, cutoff; pos=1, neg=1)
    vmin = [-1 -1 -1] 
    vmax = [1 1 1]
    invCell = inv(cell)

    for d = 1:3 
       icell = norm(invCell[d,:])
       l = floor(cutoff * icell)
       vmin[d] = -(neg + l)
       vmax[d] = pos + l 
    end 
    return get_latvec_ids(vmin, vmax)
end

function get_tbcells(atoms, invcell, mesh, ia)
    tbcells = zeros(Int64, length(atoms), 3)
    for ja = 1:length(atoms)
        Rij = get_R(atoms, ia, ja)
        tbcells[ja,:] = convert(Array{Int64,1}, round.(invcell * Rij))
    end
    return tbcells
end

function wrap_shift(shift, mesh)
    m = mesh .÷ 2
    for a = 1:length(shift)
        if shift[a] < - m[a]
            shift[a] += mesh[a]
        elseif a > m[a]
            shift[a] -= mesh[a]
        end
    end
    return shift
end

function get_neighbours(atoms, tbcells, cutoff)

    nNeighs = [ 0 for i=1:length(atoms)]
    iNeighs = [[] for i=1:length(atoms)]
    img2CentCell = []
    iCellVec = []
    coords = []
    species = []
    cell = atoms.cell
    sqcut = cutoff^2

    nAllAtom = 0
    # Outer loop: All atoms in all translation cells
    for c = 1:length(tbcells)
        oldAtom_j = 0
        for j = 1:length(atoms)
            vec = tbcells[c] * cell
            neigh_pos = atoms.X[j][:] + vec[:]
            n = 0
            # Inner loop: All atoms in central cell
            for i = 1:length(atoms)
                dist = norm(neigh_pos .- atoms.X[i])
                if (dist > sqcut)
                    continue
                end
                if j != oldAtom_j
                   nAllAtom += 1
                   push!(coords, neigh_pos)
                   push!(img2CentCell, j)
                   push!(species, atoms.Z[j])
                   push!(iCellVec, c)
                   oldAtom_j = j
                end

                # Do not add self
                if (c == 1) && (i == j)
                   continue
                end

                nNeighs[i] += 1
                push!(iNeighs[i],nAllAtom)
            end
        end
    end
    return nNeighs, iNeighs, img2CentCell, iCellVec, coords, species
end

function get_sparse_indexing(atoms, nNeighs, iNeighs, img2CentCell, norbs)

    nAtoms = size(iNeighs,1)
    mNeighs = size(iNeighs,2)
    iPairs= [[] for i=1:nAtoms+1]
    all_types = [ atoms.Z[i].z for i=1:length(atoms)]
    types = unique!(sort!(all_types))
    a2t = [findfirst(isequal(atoms.Z[i].z), types) for i=1:length(atoms) ]

    ind = 1
    for i = 1:nAtoms
      isp = a2t[i]

      # Add onsite block starting index
      push!(iPairs[i],ind)
      ind += norbs[isp] * norbs[isp]

      # Add off-site cell blocks starting indexes for each neighbour
      for inei = 1:nNeighs[i]
        jsp = a2t[img2CentCell[iNeighs[i][inei]]]
        push!(iPairs[i],ind)
        ind += norbs[isp] * norbs[jsp]
      end

    end
    return iPairs
end

end # End of module
