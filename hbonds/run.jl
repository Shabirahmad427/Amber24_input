include("./hbond.jl")
sim = load_simulation("shabir.pdb", "shabir.nc")
protein_hbonding(sim)

##julia -t 16 --project=. run.jl
