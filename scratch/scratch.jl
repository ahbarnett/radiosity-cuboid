# initial attempt to template the struct, kept havaing type errors...

# mutable struct Face{T}      # doesn't need to be mutable if v changes
# 	v::Vector{SVector{3,T}}   # vertices in RH rule order about n
# 	n::SVector{3,T}       # normal
# 	sig::Matrix{SVector{3,T}}    # radiosity values on grid
# 	w::Matrix{SVector{3,T}}    # quadrature wei on grid
# 	function Face(v::T) where T        # build unit normal via first 3 pts
# 		n = cross(v[2]-v[1],v[3]-v[2])   # assume not collinear
# 		new{T}(v, n/norm(n))
# 	end
# end

