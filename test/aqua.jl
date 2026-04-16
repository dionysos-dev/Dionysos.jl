module TestAqua

using Test
using Aqua
using Dionysos

@testset "Aqua" begin
    Aqua.test_all(Dionysos)
end

end
