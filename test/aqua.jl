module TestAqua

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Aqua

@testset "Aqua" begin
    Aqua.test_all(Dionysos)
end

end
