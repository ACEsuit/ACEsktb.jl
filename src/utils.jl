module Utils

import JuLIP: Atoms
import HDF5, JSON

function get_data(filenames, cutfunc)
    data = []
    for f in filenames
        SKdata = h5read_SK(f; get_HS=true, 
                              get_atoms=true, 
                              get_metadata=true, 
                              get_energies=true)

        kk = size(SKdata[1],2)
        ijmat = Int.(SKdata[1][1:2,1:kk])
        Rijmat = SKdata[1][3:5,1:kk]
        sel_ijR = [1,3,4,5]
        ijRijmat = SKdata[1][sel_ijR,1:kk]
        Bondint = SKdata[1][6:end,1:kk]
        at = SKdata[6]

        for i in 1:size(Bondint,2)
            R0 = SVector(Rijmat[:,i]...)
            ii = ijmat[1,i]
            Renv = get_env(at, R0, ii, cutfunc)
            E = Bondint[:,i]
            push!(data,[R0,Renv,E])
        end
    end
    return data
end

function read_json(fname)
  return JSON.parse(join(readlines(fname)))
end

function h5read_SK(fname; get_HS=false, get_atoms=false, get_metadata=false, get_energies=false)
    data = []
    HDF5.h5open(fname, "r") do fd
        groupnames = []
        for obj in fd
           push!(groupnames, HDF5.name(obj))
        end
        group_name = groupnames[1]
        dataset_names = []
        for obj in fd[group_name]
           push!(dataset_names, HDF5.name(obj))
        end
        SK_HS = HDF5.read(fd, string(group_name,"/SKHS"))
        push!(data, SK_HS)
        if get_HS
            HS_data = []
            H_str = string(group_name,"/H")
            S_str = string(group_name,"/S")
            if haskey(fd, H_str)
                H = HDF5.read(fd, string(group_name,"/H"))
                push!(HS_data, H)
            end
            if haskey(fd, S_str)
                S = HDF5.read(fd, string(group_name,"/S"))
                push!(HS_data, S)
            end
            push!(data, HS_data)
        end
        if get_energies
            energy = HDF5.read(fd, string(group_name,"/energy"))
            freeenergy = HDF5.read(fd, string(group_name,"/freeenergy"))
            push!(data, [energy, freeenergy])
        end
        if get_atoms
            positions = HDF5.read(fd, string(group_name,"/positions"))
            unitcell = HDF5.read(fd, string(group_name,"/unitcell"))
            species = HDF5.read(fd, string(group_name,"/species"))
            push!(data, [unitcell, species, positions])
        end
        if get_metadata
            metadata_str = HDF5.read(fd, string(group_name,"/metadata"))
            metadata = JSON.parse(metadata_str)
            push!(data, metadata)
        end
        push!(data, JuLIP.Atoms(; X = positions, Z = species, 
                                  cell = unitcell,
                                  pbc = [true, true, true]))
    end
    return data
end

end # End of Module
