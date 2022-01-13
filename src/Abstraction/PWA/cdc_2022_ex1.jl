using Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra

const AB = Dionysos.Abstraction

lib = CDDLib.Library() #polyhedron lib
# aux functions
eye(n) = diagm(ones(n)) # I matrix
sm(M) = SMatrix{size(M,1),size(M,2)}(M)
sv(M) = SVector{size(M,1)}(M)

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

      is_controllable, K0, P0, gamma = AB._provide_P(system);

      vol_P0 = AB.ellipsoid_vol(P0,1)
      P = P0*(vol_P0/initial_vol)^(2/n_sys)
      #println(round.(P,digits=4))
      Pp = P*contraction
      c = hcat([0.0;0.0;0.0]);
      #cp =B*5
      cp =hcat([0.1;0.5;1.9]);
      if (cp'P*cp)[1] ≤ 1
            println("cp in B_s")
      end

      @time has_transition, cost, kappa = AB._has_transition(system,P,c,Pp,cp,W,L,U)
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

contraction = 1.0;

initial_vol_span = 1:2:100 #10 .^ (0:0.1:2);
Wmax_span = 0:0.01:0.05;
contraction_span = 0.9:0.01:1.2;

cost_vector = zeros(length(initial_vol_span),length(contraction_span));
sr_vector = zeros(length(initial_vol_span),length(contraction_span));

cost_dist_vector = zeros(length(initial_vol_span),length(Wmax_span));
sr_dist_vector = zeros(length(initial_vol_span),length(Wmax_span));

for i=1:length(initial_vol_span)
      for j=1:length(contraction_span)
            has_transition, cost, sr = trial(dt,Usz,Wmax,contraction_span[j],initial_vol_span[i])
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
            has_transition, cost, sr = trial(dt,Usz,Wmax_span[j],contraction,initial_vol_span[i])
            if has_transition
                  cost_dist_vector[i,j] = cost
                  sr_dist_vector[i,j] = sr
            else
                  cost_dist_vector[i,j] = Inf
                  sr_dist_vector[i,j] = Inf
            end
      end
end


using PyPlot
PyPlot.pygui(true) 

PyPlot.rc("text",usetex=true)
PyPlot.rc("text.latex",preamble="\\usepackage{amsfonts}")

fig = PyPlot.figure(tight_layout=true, figsize=(3,3))
ax = PyPlot.axes()
#ax.set_xscale("log")
CS = PyPlot.contour(initial_vol_span,contraction_span,cost_vector')
PyPlot.clabel(CS, inline=1, fontsize=10)
PyPlot.xlabel("\${\\rm vol}(\\mathbb{B}_s)\$", fontsize=14)
PyPlot.ylabel("\$\\eta\$", fontsize=14)
PyPlot.title("\$\\rho(A_{\\rm cl})\$", fontsize=14)
plt.savefig("ex1_cost.eps", format="eps")




fig = PyPlot.figure(tight_layout=true, figsize=(3,3))
ax = PyPlot.axes()
#ax.set_xscale("log")
CS = PyPlot.contour(initial_vol_span,contraction_span,sr_vector')
PyPlot.clabel(CS, inline=1, fontsize=10)
PyPlot.xlabel("vol")
PyPlot.ylabel("contraction")
PyPlot.title("specrad")

plt.savefig("ex1_sr.eps", format="eps")



fig = PyPlot.figure(tight_layout=true, figsize=(3,3))
ax = PyPlot.axes()
#ax.set_xscale("log")
CS = PyPlot.contour(initial_vol_span,Wmax_span,cost_dist_vector')
PyPlot.clabel(CS, inline=1, fontsize=10)
PyPlot.xlabel("vol")
PyPlot.ylabel("ommax")
PyPlot.title("cost")
plt.savefig("ex1_cost_omega.eps", format="eps")





fig = PyPlot.figure(tight_layout=true, figsize=(3,3))
ax = PyPlot.axes()
#ax.set_xscale("log")
CS = PyPlot.contour(initial_vol_span,Wmax_span,sr_dist_vector')
PyPlot.clabel(CS, inline=1, fontsize=10)
PyPlot.xlabel("\${\\rm vol}(\\mathbb{B}_s)\$", fontsize=14)
PyPlot.ylabel("\$\\omega_{\\max}\$", fontsize=14)
PyPlot.title("\$\\rho(A_{\\rm cl})\$", fontsize=14)
plt.savefig("ex1_sr_omega.eps", format="eps")





# ax = PyPlot.subplot(2,1,1)
# ax.semilogx(initial_vol_span,cost_vector)
# ax = PyPlot.subplot(2,1,2)
# ax.semilogx(initial_vol_span,sr_vector)
