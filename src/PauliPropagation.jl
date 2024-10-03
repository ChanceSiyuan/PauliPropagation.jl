module PauliPropagation

include("datatypes.jl")

include("gates.jl")
export
    Gate,
    PauliGate,
    FastPauliGate,
    StaticGate,
    tofastgates

include("circuits.jl")
export
    bricklayertopology,
    get2dstaircasetopology,
    hardwareefficientcircuit,
    tfitrottercircuit,
    heisenbergtrottercircuit,
    su4ansatz,
    qcnnansatz,
    appendSU4!

include("./PauliAlgebra/PauliAlgebra.jl")
export
    inttosymbol,
    symboltoint,
    getelement,
    setelement!,
    show,
    containsXorY,
    containsYorZ

include("apply.jl")
export apply  # What should I export here?

include("truncations.jl")

include("Propagation/Propagation.jl")
export mergingbfs

include("stateoverlap.jl")
export evalagainstzero, evalagainsinitial

end
