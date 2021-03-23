module TBhelpers

using ACEtb.SlaterKoster

export set_SK_orbitals, find_row_in_matrix, get_latvec_ids, get_R, get_tbcells, wrap_shift

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

end # End of module
