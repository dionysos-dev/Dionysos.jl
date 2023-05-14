using Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra
using Plots
using SDPA, JuMP

# Example to illustrate the cost of a transition between two ellispoids as a function of meta-parameters
# such as the volume of the initial ellispoid and the contraction factor


lib = CDDLib.Library() #polyhedron lib
# aux functions
eye(n) = diagm(ones(n)) # I matrix
sm(M) = SMatrix{size(M,1),size(M,2)}(M)
sv(M) = SVector{size(M,1)}(M)

optimizer = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)

function trial(dt,Usz,Wmax,contraction,initial_vol)
      # Define system
      Ac = sm([0.0  1.0  0.0;
               0.0  0.0  1.0;
               1.0 -1.0 -1.0]);

      Bc = sm([0.0; 0.0; 1.0]);

      gc = sv(zeros(3,1))

      n_sys = 3
      n_u = 1; 

      Ze = zeros(n_u+1,n_sys+n_u+1)

      M_aux = exp(vcat([Ac Bc gc],Ze)*dt)


      A = M_aux[1:n_sys,1:n_sys];

      B = M_aux[1:n_sys,n_sys.+(1:n_u)]
      g = M_aux[1:n_sys,n_sys+n_u+1];

      Uaux = diagm(1:n_u)
      U = [(Uaux.==i)./Usz for i in 1:n_u];




      W = Wmax*hcat(collect.(vec(collect(Iterators.product(eachrow(repeat([-1 1],n_sys))...))))...); # polytope of disturbances

      L = [eye(n_sys+n_u) zeros(n_sys+n_u,1)]*dt;




      repX = intersect(HalfSpace(SVector{n_sys}(zeros(n_sys,1)), 0))
      pX = polyhedron(repX, lib);

      repU = HalfSpace(SVector{n_u}(-[1.0]), Usz) ∩ HalfSpace(SVector{n_u}([1.0]), Usz) 
      pU = polyhedron(repU, lib);

      system = ConstrainedAffineControlDiscreteSystem(A, B, g, pX, pU) 

      is_controllable, K0, P0, gamma = Dionysos.Symbolic._provide_P(system, optimizer);

      vol_P0 = Dionysos.Symbolic.ellipsoid_vol(P0,1)
      P = P0*(vol_P0/initial_vol)^(2/n_sys)
      #println(round.(P,digits=4))
      Pp = P*contraction
      c = hcat([0.0;0.0;0.0]);
      #cp =B*5
      cp =hcat([0.1;0.5;1.9]);
      if (cp'P*cp)[1] ≤ 1
            println("cp in B_s")
      end

      has_transition, cost, kappa = Dionysos.Symbolic._has_transition(system,P,c,Pp,cp,W,L,U,optimizer)
      K = kappa[:,1:n_sys];
      ell = kappa[:,n_sys+1];
      sr = max(abs.(eigen(A+B*K).values)...);
      println("Has transition: $(has_transition)")
      if has_transition
            println("K:\t $(K)\nell:\t $(ell)")
            println("cost:\t $(cost)")
            println("s.r.:\t $(sr)")
      end
      return has_transition, cost, sr
end


#to vary 
# - initial volume
# - dt
# - contraction 

dt = 0.5

Usz = 20 # upper limit on |u|

Wmax = .01
initial_vol = 10

contraction = 0.8; #1.0

initial_vol_span = 1:2:100#100 #10 .^ (0:0.1:2);
Wmax_span = 0:0.01:0.05;
contraction_span = 0.9:0.01:1.0#1.2;

cost_vector = zeros(length(initial_vol_span),length(contraction_span));
sr_vector = zeros(length(initial_vol_span),length(contraction_span));

cost_dist_vector = zeros(length(initial_vol_span),length(Wmax_span));
sr_dist_vector = zeros(length(initial_vol_span),length(Wmax_span));

for i=1:length(initial_vol_span)
      for j=1:length(contraction_span)
            has_transition, cost, sr = trial(dt, Usz, Wmax, contraction_span[j], initial_vol_span[i])
            if has_transition
                  cost_vector[i,j] = cost
                  sr_vector[i,j] = sr
            else
                  cost_vector[i,j] = Inf
                  sr_vector[i,j] = Inf
            end
      end
end

for i=1:length(initial_vol_span)
      for j=1:length(Wmax_span)
            has_transition, cost, sr = trial(dt, Usz, Wmax_span[j], contraction, initial_vol_span[i])
            if has_transition
                  cost_dist_vector[i,j] = cost
                  sr_dist_vector[i,j] = sr
            else
                  cost_dist_vector[i,j] = Inf
                  sr_dist_vector[i,j] = Inf
            end
      end
end


levels = [10,12,14,16,18,20,22]
c = contour(initial_vol_span, contraction_span, cost_vector', levels=levels, color=:viridis, clabels=true, cbar=true, lw=3, 
xtickfontsize=10, ytickfontsize=10, guidefontsize=18, titlefontsize=19)
xlabel!("vol \$(\\mathbb{B}_s)\$")
ylabel!("\$\\eta\$")
title!("\$\\widetilde{\\mathcal{J}}\$")
display(c)
# savefig("ex1_cost.png")


levels = 6
c = contour(initial_vol_span, Wmax_span, sr_dist_vector', levels=levels, color=:viridis, clabels=true, cbar=true, lw=3, 
xtickfontsize=10, ytickfontsize=10, guidefontsize=18, titlefontsize=19)
xlabel!("\${\\rm vol}(\\mathbb{B}_s)\$")
ylabel!("\$\\omega_{\\max}\$")
title!("\$\\rho(A_{\\rm cl})\$")
display(c)
# savefig("ex1_sr_omega.png")


levels = 6
c = contour(initial_vol_span, contraction_span, sr_vector', levels=levels, color=:viridis, clabels=true, cbar=true, lw=3, 
xtickfontsize=10, ytickfontsize=10, guidefontsize=18, titlefontsize=19)
xlabel!("vol \$(\\mathbb{B}_s)\$")
ylabel!("\$\\eta\$")
title!("\$\\rho(A_{\\rm cl})\$")
display(c)
# savefig("ex1_sr.png")


levels = 6
c = contour(initial_vol_span, Wmax_span, cost_dist_vector', levels=levels, color=:viridis, clabels=true, cbar=true, lw=3, 
xtickfontsize=10, ytickfontsize=10, guidefontsize=18, titlefontsize=19)
xlabel!("vol \$(\\mathbb{B}_s)\$")
ylabel!("\$\\omega_{max}\$")
title!("\$\\widetilde{\\mathcal{J}}\$")
display(c)
# savefig("ex1_cost_omega.png")



