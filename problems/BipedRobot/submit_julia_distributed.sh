cd problems/BipedRobot

export DIONYSOS_ROBOT_MODEL=6D
export DIONYSOS_NWORKERS=8
export DIONYSOS_NPARTS=8
export DIONYSOS_PARTITION_STRATEGY=contiguous
export DIONYSOS_SIMPLIFY=3.0

julia --project=. execution/julia_distributed/run_abstraction.jl