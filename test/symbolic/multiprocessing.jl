module TestDistributedAbstraction

using Test
using Dionysos
using MathematicalSystems
using StaticArrays
using Distributed

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic

@everywhere begin
    struct ToyAddDynamics <: Function end
    (f::ToyAddDynamics)(x, u) = x + u
end

function build_toy_abstraction()
    X = UT.HyperRectangle(SVector(-1.0), SVector(1.0))
    U = UT.HyperRectangle(SVector(-1.0), SVector(1.0))

    sys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        Main.ToyAddDynamics(),
        1,
        1,
        X,
        U,
    )
    approx = ST.DiscreteTimeCenteredSimulation(sys)

    xgrid = MP.GridFree(SVector(0.0), SVector(1.0))
    ugrid = MP.GridFree(SVector(0.0), SVector(1.0))

    Xmap = MP.ExplicitGridMapping{1, Float64}(xgrid, X, MP.CENTER)
    Umap = MP.ExplicitGridMapping{1, Float64}(ugrid, U, MP.CENTER)

    Xset = MP.MappingSet{1}()
    Rset = MP.MappingSet{1}()
    Uset = MP.MappingSet{1}()

    sym = SY.SymbolicModelList(Xmap, Umap; Xset = Xset, Rset = Rset, Uset = Uset)

    return (; X, U, sys, approx, xgrid, ugrid, Xmap, Umap, Xset, Rset, Uset, sym)
end

@testset "Distributed abstraction helpers" begin
    toy = build_toy_abstraction()

    Xmap = toy.Xmap
    Umap = toy.Umap
    Xset = toy.Xset
    Rset = toy.Rset
    Uset = toy.Uset
    approx = toy.approx
    sym = toy.sym

    @testset "LocalGridBasedSymbolicModel forwarding" begin
        parts = SY.partition_source_state_ids(sym, 2; strategy = :roundrobin)

        Xset_local = MP.stateset_from_states(Xmap, parts[1])
        local_sym = SY.LocalGridBasedSymbolicModel(sym, Xset_local)

        @test SY.get_state_mapping(local_sym) === SY.get_state_mapping(sym)
        @test SY.get_input_mapping(local_sym) === SY.get_input_mapping(sym)
        @test SY.get_retained_set(local_sym) === SY.get_retained_set(sym)
        @test SY.get_input_set(local_sym) === SY.get_input_set(sym)

        @test collect(SY.enum_source_states(local_sym)) == parts[1]

        q = first(parts[1])
        u = first(collect(SY.enum_inputs(local_sym)))

        @test SY.get_concrete_state(local_sym, q) == SY.get_concrete_state(sym, q)
        @test SY.get_concrete_input(local_sym, u) == SY.get_concrete_input(sym, u)
    end

    @testset "partition_source_state_ids" begin
        all_states = collect(SY.enum_source_states(sym))

        rr = SY.partition_source_state_ids(sym, 2; strategy = :roundrobin)
        @test sort(reduce(vcat, rr)) == sort(all_states)

        contig = SY.partition_source_state_ids(sym, 2; strategy = :contiguous)
        @test sort(reduce(vcat, contig)) == sort(all_states)

        onepart = SY.partition_source_state_ids(sym, 1; strategy = :roundrobin)
        @test length(onepart) == 1
        @test sort(onepart[1]) == sort(all_states)

        @test_throws ErrorException SY.partition_source_state_ids(sym, 0)
        @test_throws ErrorException SY.partition_source_state_ids(
            sym,
            2;
            strategy = :unknown,
        )
    end

    @testset "_run_local_partition_ids" begin
        parts = SY.partition_source_state_ids(sym, 2; strategy = :roundrobin)

        res = SY._run_local_partition_ids(sym, approx, parts[1]; print_level = 0)

        @test res isa SY.DistributedAbstractionResult
        @test res.n_source_states == length(parts[1])
        @test res.n_transitions == length(res.transitions)
    end

    @testset "collect_abstract_transitions_distributed local fallback" begin
        trans_ref, _ =
            SY.collect_abstract_transitions(sym, approx; print_level = 0, threaded = false)
        mono = Set(trans_ref)

        trans_rr, _ = SY.collect_abstract_transitions_distributed(
            sym,
            approx;
            procs = Int[],
            nparts = 2,
            partition_strategy = :roundrobin,
            print_level = 0,
            threaded_per_worker = false,
        )

        trans_contig, _ = SY.collect_abstract_transitions_distributed(
            sym,
            approx;
            procs = Int[],
            nparts = 2,
            partition_strategy = :contiguous,
            print_level = 0,
            threaded_per_worker = false,
        )

        @test mono == Set(trans_rr)
        @test mono == Set(trans_contig)
    end

    @testset "compute_abstract_system_distributed! local fallback" begin
        sym2 = SY.SymbolicModelList(Xmap, Umap; Xset = Xset, Rset = Rset, Uset = Uset)

        SY.compute_abstract_system_distributed!(
            sym2,
            approx;
            procs = Int[],
            nparts = 2,
            partition_strategy = :roundrobin,
            print_level = 0,
        )

        trans_dist = Set(SY.enum_transitions(sym2))

        trans_ref, _ =
            SY.collect_abstract_transitions(sym, approx; print_level = 0, threaded = false)

        @test trans_dist == Set(trans_ref)
    end

    @testset "distributed local fallback with threaded worker" begin
        trans, _ = SY.collect_abstract_transitions_distributed(
            sym,
            approx;
            procs = Int[],
            nparts = 2,
            threaded_per_worker = true,
            print_level = 0,
        )

        ref, _ =
            SY.collect_abstract_transitions(sym, approx; threaded = false, print_level = 0)

        @test Set(trans) == Set(ref)
    end

    @testset "distributed pmap smoke test" begin
        if isempty(workers())
            addprocs(2)
        end

        @everywhere using Dionysos
        @everywhere using MathematicalSystems
        @everywhere using StaticArrays

        procs = workers()[1:min(2, length(workers()))]

        trans, _ = SY.collect_abstract_transitions_distributed(
            sym,
            approx;
            procs = procs,
            nparts = 4,
            print_level = 0,
        )

        ref, _ = SY.collect_abstract_transitions(sym, approx; print_level = 0)

        @test Set(trans) == Set(ref)
    end
end

end
