
using DiscreteMarkovChains, SparseArrays, LinearAlgebra, ArnoldiMethod

mutable struct MarkovChain
    #conversion of symbol
    symmodel_from_mc #list
    mc_from_symmodel #dictionnary

    transition_matrix
    chain
    shannon_entropy
    SCC
    steady_state
    entropy
    entropy_volume
    entropy_bis
    entropy_SS
    entropy_SS_finale
end

function get_entropy_SS(mc::MarkovChain)
    return sum(mc.entropy_SS)
end
function is_Markov_matrix(matrix)
    dim = size(matrix)
    bool = true
    for i=1:dim[1]
        s = 0.0
        for j=1:dim[2]
            s += matrix[i,j]
        end
        if abs(s-1.0)>0.000001
            bool = false
        end
    end
    return bool
end

function get_entropy_volume(symmodel,s)
    V = get_volume(symmodel,s)
    ϵ = 0.001
    l = log2(V/(4*ϵ*ϵ))
    return l
end

function new_MarkovChain(symmodel_from_mc, mc_from_symmodel,transition_matrix,shannon_entropy,symmodel)
    bool = is_Markov_matrix(transition_matrix)
    if !bool
        println("error Markov chain")
    else
        println("bonne Markov chain")
    end
    chain = DiscreteMarkovChain(transition_matrix)
    communications = communication_classes(chain)
    n = length(shannon_entropy)
    steady_state = zeros(n)
    entropy = zeros(n)
    entropy_volume = zeros(n)
    entropy_bis = zeros(n)
    entropy_SS = zeros(n)
    entropy_SS_finale = zeros(n)
    for (i,class) in enumerate(communications[1])
        if communications[2][i]
            transition_matrix_scc = transition_matrix[class,class]
            chain_scc = DiscreteMarkovChain(transition_matrix_scc)
            p = stationary_distribution(chain_scc)
            for (j,s) in enumerate(class)
                steady_state[s] = p[j]
                entropy[s] = steady_state[s]*shannon_entropy[s]
                entropy_volume[s] = get_entropy_volume(symmodel,s)
                entropy_bis[s] = entropy[s] #steady_state[s] + entropy_volume[s]
                if steady_state[s]>0.0
                    entropy_SS[s] = -steady_state[s]*log2(steady_state[s])
                    entropy_SS_finale[s] = entropy_SS[s] + steady_state[s]*get_entropy_volume(symmodel,s)
                end
            end
            # else let value to zero: transcient class (non recurent)
        end
    end
    #println(findall(x->(x>0),steady_state))
    return MarkovChain(symmodel_from_mc, mc_from_symmodel,transition_matrix, chain, shannon_entropy, communications, steady_state, entropy,entropy_volume,entropy_bis,entropy_SS,entropy_SS_finale)
end


function conversion(symbols)
    symmodel_from_mc = symbols #get_cells(symmodel)
    mc_from_symmodel = Dict{Int, Int}()
    for (i,s) in enumerate(symmodel_from_mc)
        mc_from_symmodel[s] = i
    end
    return symmodel_from_mc, mc_from_symmodel
end

#inputs: list of symbols, pas necessairement 1,2,3,4,5...
# list de transitions
function build_Markov_Chain(symmodel)
    symbols = get_cells(symmodel)
    symmodel_from_mc, mc_from_symmodel = conversion(symbols)

    n = length(symbols)
    transition_matrix = zeros(n,n)
    shannon_entropy = zeros(n)
    entropy = zeros(n)
    for s in symmodel_from_mc
        idx_s = mc_from_symmodel[s]
        h = 0.0
        if s == 0
            transition_matrix[idx_s,idx_s] = 1.0
        else
            for (t,symbol,proba) in get_transitions_post(symmodel, s)
                idx_t = mc_from_symmodel[t]
                transition_matrix[idx_s,idx_t] = proba
                h = h-proba*log2(proba)
            end
        end
        shannon_entropy[idx_s] = h
    end
    return new_MarkovChain(symmodel_from_mc, mc_from_symmodel,transition_matrix,shannon_entropy,symmodel)
end


function build_Markov_Chain_Sparse(symmodel)
    symbols = get_cells(symmodel)
    symmodel_from_mc, mc_from_symmodel = conversion(symbols)

    n = length(symbols)
    shannon_entropy = zeros(n)
    entropy = zeros(n)
    I = []
    J = []
    V = []
    println("BUILD SPARSE MATRIX")
    transition_matrix = spzeros(n,n)
    transition_matrix_2 = zeros(n,n)
    println(length(symmodel_from_mc))
    for s in symmodel_from_mc
        h = 0.0
        if s == 0
            transition_matrix[0,0] = 1.0
        else
            idx_s = mc_from_symmodel[s]
            for (t,symbol,proba) in get_transitions_post(symmodel, s)
                idx_t = 0
                if t == 0
                    idx_t = 0
                else
                    idx_t = mc_from_symmodel[t]
                end
                transition_matrix[idx_s,idx_t] = proba
                transition_matrix_2[idx_s,idx_t] = proba
                h = h-proba*log2(proba)
            end
        end
    end
    return new_MarkovChain_Sparse(symmodel_from_mc, mc_from_symmodel,transition_matrix,transition_matrix_2,shannon_entropy,symmodel)
end

function is_Markov_matrix_Sparse(A)
    m, n = size(A)
    if m != n
        return false
    end
    rows = rowvals(A)
    values = nonzeros(A)
    for i = 1:m
        indexes = findall(x->(x==i),rows)
        vals = values[indexes]
        s = sum(vals)
        if abs(s-1.0)>0.000001
            return false
        end
     end
     return true
end

function compute_steady_state_Sparse(M,transition_matrix_2)
    println(size(M))
    n = size(M)[1]
    # L = [I-M;ones(1,n)]
    # println(size(L))
    # b = [zeros(n,1);1]
    # println(size(b))
    #
    #
    # x = L\b
    # println(sum(x))
    # a = norm(M * x - x)
    # println(a)
    # println(max(x))
    # return x

    M = Matrix(M)
    chain = DiscreteMarkovChain(M)
    communications = communication_classes(chain)
    steady_state = zeros(n)
    for (i,class) in enumerate(communications[1])
        if communications[2][i]
            transition_matrix_scc = M[class,class]
            chain_scc = DiscreteMarkovChain(transition_matrix_scc)
            p = stationary_distribution(chain_scc)
            for (j,s) in enumerate(class)
                steady_state[s] = p[j]
            end
        end
    end
    return steady_state
    # steady_state = zeros(n)
    # entropy = zeros(n)
    # entropy_volume = zeros(n)
    # entropy_bis = zeros(n)
    # for (i,class) in enumerate(communications[1])
    #     if communications[2][i]
    #         transition_matrix_scc = transition_matrix[class,class]
    #         chain_scc = DiscreteMarkovChain(transition_matrix_scc)
    #         p = stationary_distribution(chain_scc)
    #
    #     end
    # end
    #
    # decomp, history = partialschur(M, nev=1,tol=0.00000000000001)
    # decomp.eigenvalues
    # println(decomp.eigenvalues)
    # x = decomp.Q
    # x = x/sum(x)
    # println(sum(x))
    # a = norm(M * x - x)
    # println(a)
    return x

end

function new_MarkovChain_Sparse(symmodel_from_mc, mc_from_symmodel,transition_matrix,transition_matrix_2,shannon_entropy,symmodel)
    bool = is_Markov_matrix_Sparse(transition_matrix)
    if !bool
        println("error Markov chain")
    else
        println("bonne Markov chain")
    end
    chain = nothing
    communications = nothing
    n = length(shannon_entropy)
    println("COMPUTE STEADY-STATE")
    steady_state = compute_steady_state_Sparse(transition_matrix,transition_matrix_2)
    entropy = zeros(n)
    entropy_volume = zeros(n)
    entropy_bis = zeros(n)
    return MarkovChain(symmodel_from_mc, mc_from_symmodel,transition_matrix, chain, shannon_entropy, communications, steady_state, entropy,entropy_volume,entropy_bis)
end

function get_entropy(mc::MarkovChain, s::Int)
    i = mc.mc_from_symmodel[s]
    return mc.entropy[i]
end

function get_shannon_entropy(mc::MarkovChain, s::Int)
    i = mc.mc_from_symmodel[s]
    return mc.shannon_entropy[i]
end

function get_steady_state(mc::MarkovChain, s::Int)
    i = mc.mc_from_symmodel[s]
    return mc.steady_state[i]
end


function get_entropy_chain(mc::MarkovChain)
    h = 0.0
    for p in mc.steady_state
        if p > 0
            h = h - p*log2(p)
        end
    end
    return h
end

function get_highest_entropy(mc::MarkovChain)
    i = findmax(mc.entropy)
    return (mc.symmodel_from_mc[i], mc.entropy[i])
end
