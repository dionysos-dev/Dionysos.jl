module Test
using Test, StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

sleep(0.1) # used for good printing
println("Started test")

@testset "NestedDomain" begin
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))

    # Grid settings
    hx = SVector(6.0, 2.0)
    periodic_dims = SVector(1)
    periods = SVector(30.0)
    start = SVector(0.0)

    free_space = UT.LazySetMinus(X, UT.LazyUnionSetArray([obstacle]))
    function test_with_domain(domain)
        DO.add_set!(domain, free_space, DO.OUTER)
        # Create the NestedDomain
        Ndomain = DO.NestedDomain(domain)

        # Test initial state
        @test DO.get_levels(Ndomain) == 1
        @test !isempty(Ndomain)

        # Refinement: divide certain positions
        div = 3
        new_subpos1 = DO.cut_pos!(Ndomain, (2, 2), 1; div = div)
        new_subpos2 = DO.cut_pos!(Ndomain, (2, 3), 1; div = div)
        new_subpos3 = DO.cut_pos!(Ndomain, (4, 4), 2; div = div)

        # Test after cuts
        @test DO.get_levels(Ndomain) == 3  # We should now have 3 levels

        # Check that all cut positions at higher levels are activated
        for spos in new_subpos1
            @test DO.is_active(Ndomain, spos, 2)
        end
        for spos in new_subpos2
            @test DO.is_active(Ndomain, spos, 2)
        end
        for spos in new_subpos3
            @test DO.is_active(Ndomain, spos, 3)
        end

        # Check that the original coarse positions are deactivated
        @test !DO.is_active(Ndomain, (2, 2), 1)
        @test !DO.is_active(Ndomain, (2, 3), 1)
        @test !DO.is_active(Ndomain, (4, 4), 2)

        # Check dimensions match
        @test length(Ndomain.domains) == 3
        @test all(d -> d isa DO.GridDomainType, Ndomain.domains)

        # Basic consistency: still nonempty
        @test !isempty(Ndomain)

        # Check the number of cells increased
        ncells_before = DO.get_ncells(Ndomain)
        @test ncells_before > 0

        fig = plot(; aspect_ratio = :equal, legend = false)
        plot!(fig, Ndomain)
        @test isa(fig, Plots.Plot{Plots.GRBackend})
    end
    test_with_domain(DO.PeriodicDomainList(periodic_dims, periods, start, hx))
    test_with_domain(DO.DomainList(hx))
end

end
