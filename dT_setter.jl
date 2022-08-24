using FLoops
cell=ARGS[1]
    
ncores = ARGS[2]
path = string("/home/stack/",cell)
println(path)
ncores= parse(Int64, ncores)
@floop ThreadedEx(basesize = ncores) for _ in 1:1
    run(`python Simulation_8Types_20_12.py $path`)
end
