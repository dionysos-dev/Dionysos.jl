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

@testset "Distributed abstraction helpers" begin
    # ------------------------------------------------------------
    # Small toy system
    # x⁺ = x + u, centered simulation
    # ------------------------------------------------------------
    X = UT.HyperRectangle(SVector(-1.0), SVector(1.0))
    U = UT.HyperRectangle(SVector(-1.0), SVector(1.0))

    f = (x, u) -> x + u

    sys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(f, 1, 1, X, U)

    approx = ST.DiscreteTimeCenteredSimulation(sys)

    xgrid = MP.GridFree(SVector(0.0), SVector(1.0))
    ugrid = MP.GridFree(SVector(0.0), SVector(1.0))

    Xmap = MP.ExplicitGridMapping{1, Float64}(xgrid, X, MP.CENTER)
    Umap = MP.ExplicitGridMapping{1, Float64}(ugrid, U, MP.CENTER)

    Xset = MP.MappingSet{1}()
    Rset = MP.MappingSet{1}()
    Uset = MP.MappingSet{1}()

    sym = SY.SymbolicModelList(Xmap, Umap; Xset = Xset, Rset = Rset, Uset = Uset)

    @testset "LocalGridBasedSymbolicModel forwarding" begin
        parts = SY.partition_source_states(sym, 2; strategy = :roundrobin)
        local_sym = SY.LocalGridBasedSymbolicModel(sym, parts[1])

        @test SY.get_state_mapping(local_sym) === SY.get_state_mapping(sym)
        @test SY.get_input_mapping(local_sym) === SY.get_input_mapping(sym)
        @test SY.get_retained_domain(local_sym) === SY.get_retained_domain(sym)
        @test SY.get_input_domain(local_sym) === SY.get_input_domain(sym)

        @test collect(SY.enum_source_states(local_sym)) ==
              collect(MP.enum_states(parts[1], Xmap))

        q = first(collect(SY.enum_source_states(local_sym)))
        u = first(collect(SY.enum_inputs(local_sym)))

        @test SY.get_concrete_state(local_sym, q) == SY.get_concrete_state(sym, q)
        @test SY.get_concrete_input(local_sym, u) == SY.get_concrete_input(sym, u)
    end

    @testset "partition_source_states" begin
        all_states = collect(SY.enum_source_states(sym))

        rr = SY.partition_source_states(sym, 2; strategy = :roundrobin)
        rr_states = reduce(vcat, [collect(MP.enum_states(S, Xmap)) for S in rr])
        @test sort(rr_states) == sort(all_states)

        contig = SY.partition_source_states(sym, 2; strategy = :contiguous)
        contig_states = reduce(vcat, [collect(MP.enum_states(S, Xmap)) for S in contig])
        @test sort(contig_states) == sort(all_states)

        onepart = SY.partition_source_states(sym, 1; strategy = :roundrobin)
        @test length(onepart) == 1
        @test sort(collect(MP.enum_states(onepart[1], Xmap))) == sort(all_states)

        @test_throws ErrorException SY.partition_source_states(sym, 0)
        @test_throws ErrorException SY.partition_source_states(sym, 2; strategy = :unknown)
    end

    @testset "_run_local_partition" begin
        parts = SY.partition_source_states(sym, 2; strategy = :roundrobin)
        local_sym = SY.LocalGridBasedSymbolicModel(sym, parts[1])

        res = SY._run_local_partition(local_sym, approx; print_level = 0)

        @test res isa SY.DistributedAbstractionResult
        @test res.n_source_states == length(collect(SY.enum_source_states(local_sym)))
        @test res.n_transitions == length(res.transitions)
    end

    @testset "collect_abstract_transitions_distributed local fallback" begin
        mono = Set(
            SY.collect_abstract_transitions(sym, approx; print_level = 0, threaded = false),
        )

        part_rr = Set(
            SY.collect_abstract_transitions_distributed(
                sym,
                approx;
                procs = Int[],
                nparts = 2,
                partition_strategy = :roundrobin,
                print_level = 0,
                threaded_per_worker = false,
            ),
        )

        part_contig = Set(
            SY.collect_abstract_transitions_distributed(
                sym,
                approx;
                procs = Int[],
                nparts = 2,
                partition_strategy = :contiguous,
                print_level = 0,
                threaded_per_worker = false,
            ),
        )

        @test mono == part_rr
        @test mono == part_contig
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
        trans_ref = Set(
            SY.collect_abstract_transitions(sym, approx; print_level = 0, threaded = false),
        )

        @test trans_dist == trans_ref
    end

    @testset "distributed local fallback with threaded worker" begin
        trans = SY.collect_abstract_transitions_distributed(
            sym,
            approx;
            procs = Int[],
            nparts = 2,
            threaded_per_worker = true,
            print_level = 0,
        )

        ref =
            SY.collect_abstract_transitions(sym, approx; threaded = false, print_level = 0)
        @test Set(trans) == Set(ref)
    end

    @testset "distributed pmap smoke test" begin
        if isempty(workers())
            addprocs(2)
        end
        println("Workers available for distributed test: ", workers())
        @everywhere using Dionysos

        trans = SY.collect_abstract_transitions_distributed(
            sym,
            approx;
            procs = workers()[1:min(2, length(workers()))],
            nparts = 4,
            print_level = 0,
        )

        ref = SY.collect_abstract_transitions(sym, approx; print_level = 0)
        @test Set(trans) == Set(ref)
    end
end

end
