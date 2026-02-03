# crude 3D radiosity for shadowed ground-plane with levitating convex body
# defined by flat faces.
# Barnett, afternoon & evening of 2/2/26.
using LinearAlgebra, GLMakie, StaticArrays, Printf, Base.Threads

mutable struct Face      # doesn't need to be mutable if v changes
	# non templated version
 	v  # vertices in RH rule order about n
    n    # normal
	ng   # grid size per dim
	N    # num unknowns (grid pts)
 	w    # quadrature wei as grid (matrix)
	x    # quadr nodes as matrix of svecs
	Noff   # offset for dofs in this face
	function Face(v)        # build unit normal from first 3 pts
 		n = cross(v[2]-v[1],v[3]-v[2])   # assume not collinear
 		new(v, n/norm(n))
 	end
end
f = Vector{Face}()    # the set of surfaces (flat polygons)
L = 4.0; zg = -2.0    # half-size and z-loc of of groundplane
vgnd = [SVector{3}(c) for c in eachcol([-L L L -L; -L -L L L; zg zg zg zg])]
push!(f, Face(vgnd))   # gnd face always f[1]. Special that will get shadow RHS

# 8 cuboid verts as vector of "points" (svecs), Float64
vcub = vec([SVector{3}([x,y,z]) for x in [-1.0,1.0], y in [-1.0,1.0], z in [-1.0,1.0]])
# 6 cuboid faces, each with outward normal (defined by vert ordering)
push!(f, Face(vcub[[1 3 4 2]]))   # subsets of the svec vertex lists
push!(f, Face(vcub[[5 6 8 7]]))
push!(f, Face(vcub[[2 4 8 6]]))
push!(f, Face(vcub[[3 1 5 7]]))
push!(f, Face(vcub[[4 3 7 8]]))
push!(f, Face(vcub[[1 2 6 5]]))
#not helpful due to needing CCW ordering...
#push!(f, [vcub[j] for j in eachindex(vcub) if vcub[j][dir] > 0.0]
nf = length(f)

function plot(f::Face)  # plots into current makie fig
	v = vec(f.v)
	c = RGBf(0.7*rand(3)...)   # random somewhat dark color
	mesh!(Point3f.(v),[1 2 3; 3 1 4], color=c, shading=false)  # assumes a quad: two triangles
	meanv = sum(v)/length(v)
	arrows3d!(meanv, f.n; color=:red, lengthscale=0.5)
end

fig=Figure()
ax = Axis3(fig[1, 1], aspect = :data, viewmode=:fit, title="vertices and faces (random colors)")
# vertices
vcubp = Point3f.(vcub)  # verts for plotting
scatter!(vcubp, markersize = 20, depthsorting = true)
# label the verts
[text!(vcubp[j]; text=@sprintf("%d",j)) for j in eachindex(vcubp)]
plot.(f)   # show gnd & faces
display(fig)

function discr!(f::Face)   # use ng to apply Cartesian prod TR rule to rect
# careful: parallelogram flat face only
	ng = f.ng
    a = f.v[2]-f.v[1]; b = f.v[3]-f.v[2]    # principal vecs
	dA = norm(cross(a,b))   # area of face
	w1d = [0.5; ones(ng-2); 0.5] ./ (ng-1)   # trap rule on [0 1]; sum = 1
	f.w = dA * w1d*w1d'          # outer prod
	f.N = length(f.w)
	z = (0:ng-1)/(ng-1)           # list of fracs
	f.x = [f.v[1] + u*a + v*b for u in z, v in z]  # matrix of svecs
end

f[1].ng = 101  # discr nodes per side of gnd plane
discr!(f[1])      # applies to each face separately
Noff = 0             # track dof index offset
f[1].Noff = 0
for j=2:nf
	f[j].ng = 11   # nodes per side of other faces
	discr!(f[j])
	Noff += f[j-1].N
	f[j].Noff = Noff
end
N = sum(fj.N for fj in f)  # total unknowns

x = reduce(vcat, [fj.x[:] for fj in f])    # stack of all pts as svecs
w = reduce(vcat, [fj.w[:] for fj in f])    # stack wei
@printf("%d faces total. N=%d. total surface area = %.3g.\n",nf,N,sum(w))
ax2 = Axis3(fig[2, 1], aspect = :data, viewmode=:fit, title="nodes")
scatter!(x, markersize=5)

function showdens(f::Face, dens)  # plot density (radiosity) of a face to curr ax
	# dens can be any array shape
	@assert length(dens) == length(f.w)
	X = [p[1] for p in f.x]    # break up the svectors
	Y = [p[2] for p in f.x]
	Z = [p[3] for p in f.x]
	# reshape is to make dens a matrix shape that surface likes...
	C = [RGBf(d,d,d) for d in reshape(dens, size(X))]   # monochrome for now
	surface!(X, Y, Z, color=C, colorrange=(0,1), shading=false)   # no lighting :)
end
# handle vector of faces by view into a single dens dof vec...
showdens(f::Vector{Face}, dens) = [showdens(fj, view(dens,fj.Noff.+(1:fj.N))) for fj in f]

figd=Figure(backgroundcolor = RGBAf(0, 0, 0, 1))    # dens plot as brightness on [0,1] scale
axd = Axis3(figd[1, 1], aspect = :data, viewmode=:fit, title="density")
#https://stackoverflow.com/questions/78868362/setting-the-background-black-in-glmakie-julia
#axd.backgroundcolor = "black"  # not needed if do it at fig level
#n = f[1].ng; showdens(f[1], rand(n,n))  # test single face
showdens(f, (1:N)/N)   # test all faces w/ graded dens
#showdens(f,rand(N))
#showdens(f[3], rand(f[3].N))   # top face missing?
display(figd)

dinc = SVector{3}([-0.4,0.7,-1.0])  # incident light (parallel from infty)
dinc /= norm(dinc)   # unit

function inpoly2(x, v)   # x the pts (svec2), v is list of svec2
	#https://wrfranklin.org/Research/Short_Notes/pnpoly.html
	c = false
	n = length(v)
	for i=1:n
		j = (i==1) ? n : i-1
		if ((v[i][2]>x[2]) != (v[j][2]>x[2])) &&
			(x[1] < (v[j][1]-v[i][1]) * (x[2]-v[i][2])/(v[j][2]-v[i][2]) + v[i][1])
			c = !c
		end
	end
	c
end
v2test = vec([SVector{2}(v[ii]) for v in f[2].v])   # test bot face in 2D
@assert inpoly2(SVector{2}(0.0,0.0), v2test) 
@assert ~inpoly2(SVector{2}(2.0,0.0), v2test) 

function hits(f::Face, b,r)   # ray:  b = base pt, r = direc  (unnormalized svecs)
	# ray x(t) = b + t*r
	d = f.v[1]-b   # displacement to any pt in plane
	n = f.n
	t = dot(d,n) / dot(r,n)   # solve for ray param
	xh = b + t*r   # where hit
	# proj into 2d by dropping the coord closest to n
	_,imax = findmax(abs.(n))        # now complement it
	ii = fill(true,3); ii[imax] = 0    # now bool array giving surviving coords
	xh2 = SVector{2}(xh[ii])
	v2 = [SVector{2}(v[ii]) for v in f.v]   # convert to 2d
	#show(v2)
	inpoly2(xh2, v2)
end

#dens1 = [Float64(~hits(f[2], x, dinc)) for x in f[1].x]   # try all pts x on gnd face

gndshad = fill(1.0, f[1].N)       # compute gnd plane shadow fac
for j=2:nf
	for i=1:length(gndshad)
		gndshad[i] *= ~hits(f[j], f[1].x[i], dinc)
	end
end
#dens = rand(N)
rhs = zeros(N)
rhs[1:f[1].N] .= -dot(f[1].n,dinc) * gndshad   # all same cos factor
# now fill RHS with shadowing and cosine factors (switch off cos<0)...
for j=2:nf       # RHS for remaining faces
	cosfac = max(-dot(f[j].n,dinc), 0.0)
	rhs[f[j].Noff .+ (1:f[j].N)] .= cosfac
end
showdens(f, rhs)
display(figd)

# now fill nontrivial A block (>10 sec)...
Ng = f[1].N; Nc = N - Ng   # gnd + cuboid dofs
function fillKgc(f,Ng,Nc)     # fill dense gndface-from-cuboid block of K
	# (making a func allows better compilation)
	K = zeros(Ng,Nc)      # upper right blocks of A
	col = 1    # index in Agc
	for k=2:length(f)       # faces
		for j=1:f[k].N        # source index within kth face
			for i=1:Ng         # row (targ)
				dx = f[1].x[i] - f[k].x[j]  # displ = x-y  (target-from-src)
				nydd = dot(dx, f[k].n)   # source normal,  n_y . (x-y)
				# object above gndplane -> no need to test sign of targ normal
				K[i,col] = nydd<0.0 ? 0.0 : -nydd*dot(dx,f[1].n)/(pi*norm(dx)^4)
			end
			col+=1
		end
	end
	K
end
Kgc = fillKgc(f,Ng,Nc)
Kcg = Kgc'   # do this to kernel (before quad wei)
# applies off-diag blocks only (to quad-wei dens)...
applyA(dens) = [Kgc*view(w.*dens,Ng+1:N); Kcg*view(w.*dens,1:Ng)] # w=quadr wei
# sys is now (I-A).dens = rhs

#axr = Axis3(figd[2, 1], aspect = :data, viewmode=:fit, title="radiosity soln")
#axr.backgroundcolor = "black"   # didn't help remove black axes in axr
dens = copy(rhs)   # init Picard iter: 0th term in Neumann series
for k=1:10        # Picard iter
	@printf("total emitted power = %10.6g\t\tPicard iter %d...\n",sum(w.*dens),k)
	dens += applyA(dens)
end
showdens(f, dens)
# to do: try CG on I + K.   
display(figd)