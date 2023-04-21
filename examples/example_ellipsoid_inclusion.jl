include("../src/Dionysos.jl")
using .Dionysos
UT = Dionysos.Utils
using Plots, Colors, LinearAlgebra, LaTeXStrings


myblue = RGB(108 ./256,142 ./256,191 ./256) 
myblueN = RGB(32 ./256,103 ./256,205 ./256)
myorange = RGB(255 ./280,158 ./280,56 ./280) 
myorangeN = RGB(255 ./280,130 ./280,0.0 ./280) 

function tansformations()
    c0 = [3.0; 5.0]
    P0 = [3.0 0.0;0.0 1/4.0]

    c = [-4.0; 7.0]
    P = [0.4 -0.1;-0.1 0.5] #[2.0 0.0;0.0 5.0]

    El0 = UT.Ellipsoid(P0, c0)
    El = UT.Ellipsoid(P, c)

    ## Change of variables 1 ##
    L0 = cholesky(El0.P).L
    El0_1 = UT.transform(El0, L0',-L0'*El0.c)
    El_1 = UT.transform(El, L0',-L0'*El0.c) 
    ## Change of variables 2 ##
    specDecomp = eigen(El_1.P)
    V = specDecomp.vectors
    El0_2 = UT.transform(El0_1, V',0*El0.c) 
    El_2 = UT.transform(El_1, V',0*El0.c) 

    Elnew = UT.scale_for_inclusion_contact_point(El0_2, El_2)

    println(El0∈El)

    p = plot(aspect_ratio=:equal)
    UT.plotE!(El0, color=:orange, label="El0")
    UT.plotE!(El, color=:blue, label="El")

    UT.plotE!(El0_1, color=:black, label="El0_1", opacity=0.5)
    UT.plotE!(El_1, color=:yellow, label="El_1")

    UT.plotE!(Elnew, color=:pink, label="Elnew")
    UT.plotE!(El0_2, color=:green, label="El0_2")
    UT.plotE!(El_2, color=:red, label="El_2")

    display(p)
end

function filter_by_interval(x, fx, interval)
    # Find the indices where fx is within the interval
    idx = findall(interval[1] .<= fx .<= interval[2])
    # Filter x and fx by the selected indices
    return x[idx], fx[idx]
end

function set_to_zero_if_close(vec::Vector, eps::Real)
    vec[abs.(vec) .< eps] .= 0.0
    return vec
end

function plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h; eps=1e-10)
    L = cholesky(El0.P).L
    P = L\El.P/L';
    c = L'*(El.c -El0.c)
    specDecomp = eigen(P)
    vals = specDecomp.values
    ct = specDecomp.vectors'*c
    g(β) = -β + sum((β*vals./(1 .- β*vals)).*(ct.^2))

    lb = 1/min(vals...)
    ub = 1-norm(ct)^2
    println("ct: ", ct)
    gstar, βstar = UT.get_ℓ_ast_inclusion(El, El0)
    println("lb : " , lb)
    println("ub : ", ub)
    println("βstar : " , βstar)
    println("ℓstar : ", -gstar)

    p = plot(xlabel="β", ylabel="ℓ(β)", xtickfontsize=8, ytickfontsize=8, xguidefontsize=15, yguidefontsize=15, linewidth=1, title="")
    xticks!([0.0, intervalx[2]])
    yticks!([intervalfx[1],intervalfx[2]])

    inv_vals = 1 ./ vals
    tab = vcat(inv_vals, [intervalx[1], intervalx[2]])
    bounds = sort(tab)
    for i in 1:length(bounds)-1
        current_interval = [bounds[i], bounds[i+1]]
        I = current_interval[1]+ϵ:h:current_interval[2]-ϵ
        I, gI = filter_by_interval(I, g.(I), intervalfx)
        plot!(I, gI, linewidth=1, color=:blue, label="")
    end
    hline!([-1], linewidth=2, color=:red, label="")
    inv_vals = 1 ./ vals
    xlimits = xlims()
    ylimits = ylims()
    unique!(inv_vals)
    vline!(inv_vals, color=:black, linewidth=2, label="")
    for (i,x) in enumerate(inv_vals)
        annotate!(x, ylimits[1]-0.08*(ylimits[2]-ylimits[1]), text(latexstring("\$\\frac{1}{\\lambda_{$(i)}}\$")))
    end
    annotate!(xlimits[1] + 0.04*(xlimits[2]-xlimits[1]), -1+0.04*(ylimits[2]-ylimits[1]), text("-1", :right, color=:red, 14))

    hline!([-gstar], linestyle=:dash, linewidth=2, color=:green, label="")
    annotate!(xlimits[1]- 0.01*(xlimits[2]-xlimits[1]), -gstar, text("ℓ*", :right, color=:green, 14))
    vline!([βstar], linestyle=:dash, linewidth=2, color=:green, label="")
    annotate!(βstar, ylimits[1]- 0.05*(ylimits[2]-ylimits[1]), text("β*", :upper, color=:green, 14))
    scatter!([βstar], [-gstar], markersize=4, color=:green, label="")

    if ub>=0.0
        vline!([ub], linewidth=2, color=:black, label="")
        annotate!(ub, ylimits[1]- 0.05*(ylimits[2]-ylimits[1]), text("ub", :upper, color=:black, 14))
    end
    display(p)
end


function fig1()
    a = 1.0
    c0 = [1.6+a; 1.4+a]
    P0 = [0.4 -0.1;
         -0.1 0.5]
    
    c = [1.5; 1.5]
    P = [4.0 0.5;       
         0.5 6.0]
    ############ Plot parameters ############
    intervalx = [-0.02, 0.4]
    intervalfx = [-3, 1.0]
    ϵ = 0.0001
    h = 0.0001
    ########################################
    El0 = UT.Ellipsoid(P0, c0)
    El = UT.Ellipsoid(P, c)
    #########################################
    p = plot(aspect_ratio=:equal,legend=false)
    UT.plotE!(El0, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    annotate!(3.0, 3.0, text(latexstring("\$\\mathcal{E}_0\$"), :upper, color=myorangeN, 30))
    annotate!(1.5, 1.5, text(latexstring("\$\\mathcal{E}\$"), :upper, color=myblueN, 30))
    display(p)
    #########################################
    plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h)
end

function fig2()
    a = 0.89
    c0 = [1.6+a; 1.4+a]
    P0 = [0.4 -0.1;
         -0.1 0.5]
    
    c = [1.5; 1.5]
    P = [4.0 0.5;       
         0.5 6.0]
    ############ Plot parameters ############
    intervalx = [-0.02, 0.5]
    intervalfx = [-3, 1.0]
    ϵ = 0.0001
    h = 0.0001
    ########################################
    El0 = UT.Ellipsoid(P0, c0)
    El = UT.Ellipsoid(P, c)
    #########################################
    p = plot(aspect_ratio=:equal,legend=false)
    UT.plotE!(El0, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    annotate!(3.0, 3.0, text(latexstring("\$\\mathcal{E}_0\$"), :upper, color=myorangeN, 30))
    annotate!(1.5, 1.5, text(latexstring("\$\\mathcal{E}\$"), :upper, color=myblueN, 30))
    display(p)
    #########################################
    plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h)
end

function fig3()
    a = 0.6
    c0 = [1.6+a; 1.4+a]
    P0 = [0.4 -0.1;
         -0.1 0.5]
    
    c = [1.5; 1.5]
    P = [4.0 0.5;       
         0.5 6.0]
    ############ Plot parameters ############
    intervalx = [-0.02, 0.8]
    intervalfx = [-3, 1.0]
    ϵ = 0.0001
    h = 0.0001
    ########################################
    El0 = UT.Ellipsoid(P0, c0)
    El = UT.Ellipsoid(P, c)
    #########################################
    p = plot(aspect_ratio=:equal,legend=false)
    UT.plotE!(El0, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    annotate!(3.0, 3.0, text(latexstring("\$\\mathcal{E}_0\$"), :upper, color=myorangeN, 30))
    annotate!(1.5, 1.5, text(latexstring("\$\\mathcal{E}\$"), :upper, color=myblueN, 30))
    display(p)
    #########################################
    plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h)
end

function fig4()
    a = 0.6
    c0 = [1.6+a; 1.4+a]
    P0 = [0.4 -0.1;
         -0.1 0.5]
    
    c = [1.5; 1.5]
    P = [4.0 0.5;       
         0.5 6.0]
    ############ Plot parameters ############
    intervalx = [-0.02, 0.8]
    intervalfx = [-3, 1.0]
    ϵ = 0.0001
    h = 0.0001
    ########################################
    El0 = UT.Ellipsoid(P0, c0)
    El = UT.Ellipsoid(P, c)
    #########################################
    p = plot(aspect_ratio=:equal,legend=false)
    UT.plotE!(El0, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    annotate!(3.0, 3.0, text(latexstring("\$\\mathcal{E}_0\$"), :upper, color=myorangeN, 30))
    annotate!(1.5, 1.5, text(latexstring("\$\\mathcal{E}\$"), :upper, color=myblueN, 30))
    display(p)
    #########################################
    plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h)
end

function particularCase()
    c3 = [0.0; 1.5]
    P3 = [1.0 0.0;       
          0.0 2.0]
    El3 = UT.Ellipsoid(P3, c3)
    #### change of variables ####
    θ = π/3  # 30 degrés en radians
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    c2 = R*c3
    P2 = R*El3.P*R'
    El2 = UT.Ellipsoid(P2, c2)
    #### change of variables ####
    a = 0.6
    c0 = [1.6+a; 1.4+a]
    θ0 = π/4  # 30 degrés en radians
    R0 = [cos(θ0) -sin(θ0); sin(θ0) cos(θ0)]
    #P0 = [0.4 -0.1;
    #     -0.1 0.5] # c'est une matrice de réflexion}
    P0 = R0*[1.0 0.0;0.0 3.0]*R0'
    El0 = UT.Ellipsoid(P0, c0)
    El02 = UT.Ellipsoid([1.0 0.0;0.0 1.0], [0.0;0.0])
    L0 = cholesky(El0.P).L
    P1 = L0*P2*L0'
    c1 = inv(L0)'*El2.c + El0.c
    El1 = UT.Ellipsoid(P1, c1)
    #########################################
    myblue = RGB(108 ./256,142 ./256,191 ./256) 
    myblueN = RGB(32 ./256,103 ./256,205 ./256)
    myorange = RGB(255 ./280,158 ./280,56 ./280) 
    myorangeN = RGB(255 ./280,130 ./280,0.0 ./280) 

    p = plot(aspect_ratio=:equal,legend=false)
    UT.plotE!(El0, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El1, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    UT.plotE!(El02, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El2, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    UT.plotE!(El3, color=myblue; opacity=0.8, lw=4,lc=myblueN)

    UT.plot_arrow!(El02.c,El3.c;dims=[1,2],color=:black,markeralpha=1.0)
    UT.plot_arrow!(El02.c,El2.c;dims=[1,2],color=:black,markeralpha=1.0)
    UT.plot_arrow!(El0.c,El1.c;dims=[1,2],color=:black,markeralpha=1.0)

    annotate!(El0.c[1], El0.c[2], text(latexstring("\$\\mathcal{E}_0\$"), :upper, color=myorangeN, 25))
    annotate!(El1.c[1], El1.c[2], text(latexstring("\$\\mathcal{E}\$"), :upper, color=myblueN, 25))
    annotate!(El02.c[1], El02.c[2], text(latexstring("\$\\tilde{\\mathcal{E}}_0\$"), :upper, color=myorangeN, 20))
    annotate!(El2.c[1], El2.c[2], text(latexstring("\$\\tilde{\\mathcal{E}}\$"), :upper, color=myblueN, 20))
    annotate!(El3.c[1], El3.c[2], text(latexstring("\$\\tilde{\\tilde{\\mathcal{E}}}\$"), :upper, color=myblueN, 15))

    UT.plotAxis!(El0)
    UT.plotAxis!(El1)
    UT.plotAxis!(El2)
    UT.plotAxis!(El3)

    println(eigen(El1.P).vectors)
    println(eigen(El0.P).vectors)
    display(p)
   # return El1, El0
end

# pour ne pas avoir d'asymptote, il faut que l'axe avec lequel c-c0 est ortogonal soit celui avec le plus grand ratio
# la composante du vecteur c nulle sera celle de l'axe avec lequel c-c0 est orthogonale
function particularCase(c0, V0, D0, cx, i, V, D)
    ########################################
    cy = ((V[1,i]*c0[1]+V[2,i]*c0[2])-V[1,i]*cx)/V[2,i] # choisir la composante y de c pour que l'axe i de Ε soit orthogonal avec le vecteur c-c0
    c = [cx; cy]
    El0 = UT.Ellipsoid(V0*D0*V0', c0)
    El = UT.Ellipsoid(V*D*V', c)
    #########################################
    ## Change of variables 1 ##
    L0 = V0*sqrt.(D0)# cholesky(El0.P).L
    println(V0*sqrt.(D0))
    println(cholesky(El0.P).L)
    El0_1 = UT.transform(El0, L0',-L0'*El0.c)
    El_1 = UT.transform(El, L0',-L0'*El0.c) 
    ## Change of variables 2 ##
    specDecomp = eigen(El_1.P)
    V1 = V0'*V #specDecomp.vectors
    El0_2 = UT.transform(El0_1, V1',0*El0.c) 
    El_2 = UT.transform(El_1, V1',0*El0.c) 
    #########################################
    UT.plotE!(El0, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    UT.plotE!(El0_1, color=myorange; opacity=0.6, lw=4,lc=myorangeN)
    UT.plotE!(El_1, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    UT.plotE!(El_2, color=myblue; opacity=0.8, lw=4,lc=myblueN)
    UT.plot_arrow!(El0.c,El.c; markeralpha=1.0)
    UT.plot_arrow!(El0_1.c,El_1.c; markeralpha=1.0)
    UT.plotAxis!(El0;color1=:green, color2=:red)
    UT.plotAxis!(El;color1=:green, color2=:red)
    UT.plotAxis!(El_1;color1=:green, color2=:red)
    UT.plotAxis!(El_2;color1=:green, color2=:red)
    #########################################
    vals = eigen(El.P).values
    vals0 = eigen(El0.P).values
    println("ratio 1 :" ,vals0[1]/vals[1])
    println("ratio 2 :" ,vals0[2]/vals[2])
    println("Les valeurs propres: ", eigen(El_2.P).values)
    println("Le vecteur c : ", El_2.c)
    vals = eigen(El_2.P).values
    (vals[1] < vals[2] && abs.(El_2.c[1])) < 10e-10 || (vals[1] > vals[2] && abs.(El_2.c[2])) ? println("N'a pas d'asymptote") : println("A une asymptote")
    return El0, El, El0_1, El_1, El_2
end

# plot a case where the longest axis ofE is orthogonal with c-c0 but we still have an asymptote
function test5()
    c0 = [3.2; 3.0]
    θ0 = 30*π/180
    V0 = [cos(θ0) -sin(θ0); sin(θ0) cos(θ0)]
    D0 = [2.8 0.0; 0.0 10.0]

    θ = 30*π/180
    V = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    D = [1.0 0.0; 0.0 2.0]

    i = 1 # indique l'axe de Epsilon avec lequel c-c0 est orthogonal
    cx = 4.0 # composante x de c

    #########################################
    r = 20
    p = plot(aspect_ratio=:equal,legend=false)
    El0, El, El0_1, El_1, El_2 = particularCase(c0, V0, D0, cx, i, V, D)
    annotate!(El0.c[1]-0.9, El0.c[2]+0.5, text(latexstring("\$\\mathcal{E}_0\$"), color=myorangeN, r))
    annotate!(El.c[1]+1.2, El.c[2]-0.5, text(latexstring("\$\\mathcal{E}\$"), color=myblueN, r))
    annotate!(El0_1.c[1]-0.7, El0_1.c[2]-1.5, text(latexstring("\$\\tilde{\\mathcal{E}}_0\$"), color=myorangeN, r))
    annotate!(El_1.c[1]-2.1, El_1.c[2]-1.0, text(latexstring("\$\\tilde{\\mathcal{E}}\$"), color=myblueN, r))
    annotate!(El_2.c[1], El_2.c[2]-2.2, text(latexstring("\$\\bar{\\mathcal{E}}\$"), color=myblueN, r))

    annotate!(-0.25, -0.25, text(latexstring("\$0\$"), color=:black, 12))
    annotate!(El.c[1]-0.15, El.c[2]-0.3, text(latexstring("\$c\$"), color=:black, 13))
    annotate!(El0.c[1]-0.081, El0.c[2]-0.4, text(latexstring("\$c_0\$"), color=:black, 13))
    annotate!(El_1.c[1]-0.08, El_1.c[2]-0.4, text(latexstring("\$\\tilde{c}\$"), color=:black, 13))
    annotate!(El_2.c[1]-0.25, El_2.c[2]-0.3, text(latexstring("\$\\bar{c}\$"), color=:black, 13))
    display(p)
    #########################################
    intervalx = [-0.02, 40]
    intervalfx = [-200, 1.0]
    ϵ = 0.00001
    h = 0.0001
    #plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h)
end

# plot a case where we do ot have an asymptot
function test6()
    c0 = [3.2; 3.0]
    θ0 = 30*π/180
    V0 = [cos(θ0) -sin(θ0); sin(θ0) cos(θ0)]
    D0 = [6.0 0.0; 0.0 10.0]

    θ = 50*π/180
    V = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    D = [1.0 0.0; 0.0 2.0]

    i = 1 # indique l'axe de Epsilon avec lequel c-c0 est orthogonal
    cx = 4.0 # composante x de c

    #########################################
    r = 20
    p = plot(aspect_ratio=:equal,legend=false)
    El0, El, El0_1, El_1, El_2 = particularCase(c0, V0, D0, cx, i, V, D)
    annotate!(El0.c[1]-0.9, El0.c[2]+0.5, text(latexstring("\$\\mathcal{E}_0\$"), color=myorangeN, r))
    annotate!(El.c[1]+1.2, El.c[2]-0.5, text(latexstring("\$\\mathcal{E}\$"), color=myblueN, r))
    annotate!(El0_1.c[1]-0.7, El0_1.c[2]-1.5, text(latexstring("\$\\tilde{\\mathcal{E}}_0\$"), color=myorangeN, r))
    annotate!(El_1.c[1]-2.1, El_1.c[2]-1.0, text(latexstring("\$\\tilde{\\mathcal{E}}\$"), color=myblueN, r))
    annotate!(El_2.c[1], El_2.c[2]-2.2, text(latexstring("\$\\bar{\\mathcal{E}}\$"), color=myblueN, r))

    annotate!(-0.25, -0.25, text(latexstring("\$0\$"), color=:black, 12))
    annotate!(El.c[1]-0.15, El.c[2]-0.3, text(latexstring("\$c\$"), color=:black, 13))
    annotate!(El0.c[1]-0.081, El0.c[2]-0.4, text(latexstring("\$c_0\$"), color=:black, 13))
    annotate!(El_1.c[1]-0.08, El_1.c[2]-0.4, text(latexstring("\$\\tilde{c}\$"), color=:black, 13))
    annotate!(El_2.c[1]-0.25, El_2.c[2]-0.3, text(latexstring("\$\\bar{c}\$"), color=:black, 13))
    display(p)
    #########################################
    intervalx = [-0.02, 40]
    intervalfx = [-200, 1.0]
    ϵ = 0.00001
    h = 0.0001
   # plotSecularFunction(El0, El, intervalx, intervalfx, ϵ, h)
end

#tansformations()
# fig3()
test6()
# particularCase()