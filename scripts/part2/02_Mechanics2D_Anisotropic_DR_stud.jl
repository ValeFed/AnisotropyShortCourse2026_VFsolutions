using Plots, ForwardDiff, LinearAlgebra, StaticArrays

av(η) = 1/4 .* ( η[1:end-1,1:end-1] .+ η[1:end-1,2:end-0] .+ η[2:end-0,1:end-1] .+ η[2:end-0,2:end-0] )


function AnisotropicViscousRheology(ε̇_cart, η, θ, δ)

    # @TODO: Define the transformation matrix
    T     = [cosd(θ)^2 sind(θ)^2  2*cosd(θ)*sind(θ);
             sind(θ)^2  cosd(θ)^2 -2*cosd(θ)*sind(θ);
            -cosd(θ)*sind(θ) cosd(θ)*sind(θ) cosd(θ)^2 - sind(θ)^2]

    # Deviatoric strain rate in material coordinates 
    # @TODO: forward transformation
    ε̇_mat    = T * ε̇_cart

    # @TODO: Define the constitutive operator for transverse isotropic viscous flow 
    𝐃_mat    = 2η*I(3)
    𝐃_mat[3,3] = 𝐃_mat[3,3]/δ

    # @TODO: Evaluate deviatoric stress
    τ_mat = 𝐃_mat * ε̇_mat

    # Deviatoric stress in Cartesian coordinates
    # @TODO: backward transformation
    τ_cart   = T \ τ_mat
    return τ_cart, τ_mat
end

function ConstitutiveTensorVoigt(η, θ, δ)
    ε̇ = @SVector [1.0, 1.0, 1.0] 
    # @TODO: Call your function here !!!
    f(ε̇ ) = first(AnisotropicViscousRheology(ε̇ , η, θ, δ))
    𝐃 =  ForwardDiff.jacobian(f, ε̇ ) #cartesian
    return 𝐃 
end

# @TODO: Add the function you wrote in part 1 exercise 04 ...

function DeviatoricStress!(τ, ε̇, 𝐃)
    # Loop on node types (p=1 --> centroids or p=2 --> vertices) 
    for p = 1:2
        # Loop on each entry of the arrays
        for i in eachindex(𝐃[p])
            # Create a local strain rate array in Voigt notation
            e = SA[ε̇.xx[p][i]; ε̇.yy[p][i]; ε̇.xy[p][i]]
            # Apply the constitutive law
            s = 𝐃[p][i] * e
            # Unpack the stress components to the global solution arrays
            τ.xx[p][i], τ.yy[p][i], τ.xy[p][i] = s
        end
    end
end

function Mechanics2D()

    # Physics
    xmin    =-1/2
    xmax    = 1/2
    ymin    =-1/2
    ymax    = 1/2                       
    ηinc    = 1e-1
    ηmat    = 1.0
    ξmat    = 1e2
    r       = 0.1
    ε̇bg     = 1.0

    # Boundary loading type
    L_BC = SA[-ε̇bg 0.0; 
               0.0 ε̇bg]
 
    # @TODO: define transverse isotropy parameters
    δ       = 10
    θ       = 20

    # Numerics
    ncx     = 101                               # centroids in x
    ncy     = ncx                               # centroids in y
    Δx      = (xmax-xmin) / ncx                 # grid step x
    Δy      = (ymax-ymin) / ncx                 # grid step y

    # Allocate
    size_c  = (ncx+0, ncx+0)
    size_v  = (ncx+1, ncx+1)
    P       = zeros(size_c)
    η       = ( c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1) )
    ξ       = zeros(size_c)
    divV    = ( c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1) )
    V       = (  
        x  = ( i = zeros(ncx+1, ncy+2), j = zeros(ncx+2, ncy+1) ),
        y  = ( i = zeros(ncx+1, ncy+2), j = zeros(ncx+2, ncy+1) ),    
        I  = ( i = zeros(ncx+1, ncy+2), j = zeros(ncx+2, ncy+1) ), 
    )
    ε̇       = ( 
        xx = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        yy = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        zz = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        xy = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        II = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
    )
    τ       = ( 
        xx = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        yy = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        zz = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        xy = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
        II = (c = zeros(ncx+0, ncx+0), v = zeros(ncx+1, ncx+1)),
    )
    ∂Vx∂τ   = zeros(ncx+1, ncy+2) # Pseudo rate of change Vx (grid i)
    ∂Vy∂τ   = zeros(ncx+2, ncy+1) # Pseudo rate of change Vy (grid j)
    Rx      = zeros(ncx+1, ncy+2) # Momentum residual x (grid i)
    Ry      = zeros(ncx+2, ncy+1) # Momentum residual y (grid j)
    Rx0     = zeros(ncx+1, ncy+2) # Momentum residual x (grid i) old
    Ry0     = zeros(ncx+2, ncy+1) # Momentum residual y (grid j) old

    # Rheological operators for deviatoric stress-strain
    𝐃 =  ( 
        c = [@SMatrix(zeros(3,3)) for _ in 1:size_c[1], _ in 1:size_c[2]],
        v = [@SMatrix(zeros(3,3)) for _ in 1:size_v[1], _ in 1:size_v[2]],
    )

    # Pre-calculation  
    x = ( c=LinRange(xmin, xmax, ncx), v=LinRange(xmin, xmax, ncx+1), j=LinRange(xmin-Δx/2, xmax+Δx/2, ncx+2) ) 
    y = ( c=LinRange(ymin, ymax, ncy), v=LinRange(ymin, ymax, ncy+1), i=LinRange(ymin-Δy/2, ymax+Δy/2, ncy+2) ) 

    # @TODO: Define viscosity on vertices
    η.v .= ηmat
    η.v[(x.v.^2 .+ y.v'.^2) .< r^2] .= ηinc
    
    # Interpolate viscosity from vertices to centroids
    η.c .= av(η.v) 
    

    # Bulk viscosity
    ξ .= ξmat 

    # Fill in the rheological operator for each node types 
    # (p=1 --> centroids or p=2 --> vertices) 
    for p = 1:2
        for i in eachindex(𝐃[p])
            # @TODO: call the function ConstitutiveTensorVoigt()
            𝐃[p][i] = ConstitutiveTensorVoigt(η[p][i],θ,δ) #each cell has its own consitutive operator
        end
    end

    # @TODO: define a pure shear velocity Vx array on i grid and Vy array on j grid
    V.x.i  .= L_BC[1,1] * x.v # just because it is pure shear, we don't need the shear component
    V.y.j  .= L_BC[2,2] * y.v'

    # @TODO: define these scalar values of boudary velocity components
    VxW, VxE = V.x.i[1], V.x.i[end] # or L_BC[1,1]*xmin and L_BC[1,1]*xmax
    VyS, VyN = V.y.j[1], V.y.j[end] # or L_BC[2,2]*ymin and L_BC[2,2]*ymax

    # Visualise initial fields
    p1 = heatmap(x.v, y.i, V.x[1]', aspect_ratio=1, title="Vx")
    p2 = heatmap(x.j, y.v, V.y[2]', aspect_ratio=1, title="Vy")
    p3 = heatmap(x.c, y.c, η.c',    aspect_ratio=1, title="η")
    p4 = heatmap(x.c, y.c, P',      aspect_ratio=1, title="P")
    display(plot(p1, p2, p3, p4))

    # Iterative solver parameters
    niter = 10000   # Maximum number of iterations
    nout  = 200     # Check error each nout iterations
    tol   = 1e-5    # Tolerance of the solver (stop iterations when error has decreased below that level)
    CFL   = 0.95    # CFL value is 0.90 - 0.99
    cf    = 1.0     # Damping factor ~ 0.5 - 1.0
    err0  = 1.0     # Dummy initial error value

    # Maximum eigenvalue: Gershgorin circle theorem
    ηW, ηE, ηS, ηN = η.c[2:end-0,:], η.c[1:end-1,:], η.v[2:end-1,1:end-1], η.v[2:end-1,2:end-0]
    PCx   = ones(ncx+1, ncy+2)
    λxx   = @. 4 * ξmat ./ Δx .^ 2 + (ηN ./ Δy + ηS ./ Δy) ./ Δy + ηN ./ Δy .^ 2 + ηS ./ Δy .^ 2 + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx + (4 // 3) * ηE ./ Δx .^ 2 + (4 // 3) * ηW ./ Δx .^ 2
    λyx   = @. abs(ξmat ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηN ./ (Δx .* Δy)) + abs(ξmat ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηS ./ (Δx .* Δy)) + abs(ξmat ./ (Δx .* Δy) + ηN ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy)) + abs(ξmat ./ (Δx .* Δy) + ηS ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy))
    λ     = λxx .+ λyx 
    ηW, ηE, ηS, ηN = η.v[1:end-1,2:end-1], η.v[2:end,2:end-1], η.c[:,1:end-1], η.c[:,2:end-0]
    PCy   = ones(ncx+2, ncy+1)
    λmax  = maximum( λ ./ PCx[2:end-1,2:end-1] ) 
    
    # Pseudo time step is proportional to λmax (e.g., Oakley, 1995)
    Δτ    = 2 / sqrt((λmax)) * CFL 

    # Damping coefficient is proportional to λmin (e.g., Oakley, 1995)
    λmin  = 0.         
    c     = 2*sqrt.(λmin)*cf
    α     = 2 .* Δτ^2 ./ (2 .+ c.*Δτ)
    β     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)  

    # 2D mechanical solver
    for iter=1:niter

        # Keep track of previous iteration
        Rx .= Rx0
        Ry .= Ry0

        # Boundary conditions using ghost nodes
        # BC for Vx on grid i --> ∂Vx∂y
        V.x[1][:, [1, end]] .= V.x[1][:, [2, end-1]] # S/N
        # BC for Vx on grid j --> ∂Vx∂x
        V.x[2][1,   :] .= 2 * VxW .- V.x[2][2,     :] 
        V.x[2][end, :] .= 2 * VxE .- V.x[2][end-1, :] 
        # BC for Vy on grid j --> ∂Vy∂x
        V.y[2][[1, end], :] .= V.y[2][[2, end-1], :] # W/E
        # BC for Vy on grid i --> ∂Vy∂y
        V.y[1][:, 1, ] .= 2 * VyS .- V.y[1][:, 2,   ]  
        V.y[1][:, end] .= 2 * VyN .- V.y[1][:, end-1]  
        
        # Strain rate components
        divV.c  .= diff(V.x[1][:,2:end-1], dims=1) / Δx .+ diff(V.y[2][2:end-1,:], dims=2) / Δy
        divV.v  .= diff(V.x[2],            dims=1) / Δx .+ diff(V.y[1], dims=2) / Δy
        ε̇.xx.v  .= diff(V.x[2],            dims=1) / Δx .- 1/3 .* divV.v
        ε̇.yy.v  .= diff(V.y[1],            dims=2) / Δy .- 1/3 .* divV.v
        ε̇.xx.c  .= diff(V.x[1][:,2:end-1], dims=1) / Δx .- 1/3 .* divV.c
        ε̇.yy.c  .= diff(V.y[2][2:end-1,:], dims=2) / Δy .- 1/3 .* divV.c
        ε̇.xy.c  .= diff(V.x[2][2:end-1,:], dims=2) / Δy .+ diff(V.y[1][:,2:end-1], dims=1) / Δx
        ε̇.xy.v  .= diff(V.x[1],            dims=2) / Δy .+ diff(V.y[2], dims=1) / Δx

        # Pressure
        P      .= -ξ .* divV.c  
      
        # Deviatoric stress components
        DeviatoricStress!(τ, ε̇, 𝐃)

        # Discrete residual function
        Rx[2:end-1, 2:end-1] .=  (diff(τ.xx.c, dims=1) / Δx + diff(τ.xy.v[2:end-1,:], dims=2) / Δy .- diff(P, dims=1) / Δx)          
        Ry[2:end-1, 2:end-1] .=  (diff(τ.yy.c, dims=2) / Δy + diff(τ.xy.v[:,2:end-1], dims=1) / Δy .- diff(P, dims=2) / Δy)          

        # Update rate of change of solution (in "iteration" or "pseudo" time, τ)
        ∂Vx∂τ   .= Rx./PCx .+ β*∂Vx∂τ  
        ∂Vy∂τ   .= Ry./PCy .+ β*∂Vy∂τ    
  
        # Update solution
        V.x[1]     .+= α*∂Vx∂τ
        V.y[2]     .+= α*∂Vy∂τ
        
        # Interpolate to dual grid
        V.x[2][2:end-1,:] .= av(V.x[1])
        V.y[1][:,2:end-1] .= av(V.y[2])

        if mod(iter, nout)==0 || iter==1
            # Residual check 
            (iter == 1) && (err0 = norm(Rx))
            err   = norm(Rx)/err0
            @info iter, err
            err<tol && break
            isnan(err) && error("NaN error - Iterations went wrong!")
            # Rayleigh quotient (e.g., Joldes et al., 2011)
            λmin  = abs.((sum(Δτ*∂Vx∂τ.*(Rx .- Rx0)./PCx)) + (sum(Δτ*∂Vy∂τ.*(Ry .- Ry0)./PCy))) / ( sum((Δτ*∂Vx∂τ).^2) + sum((Δτ*∂Vy∂τ).^2) )
            # Dynamic evaluation of PT iteration parameters
            c = 2*sqrt.(λmin)*cf
            α = 2 .* Δτ^2 ./ (2 .+ c.*Δτ)
            β = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)
        end
        
    end

    # @TODO: calculate velocity magnitude on i grid
    V.I.i  .= sqrt.(V.x.i.^2 .+ V.y.i.^2)

    # Plane strain
    ε̇.zz.c .= @. -(ε̇.xx.c + ε̇.yy.c)
    τ.zz.c .= @. -(τ.xx.c + τ.yy.c)

    # @TODO: calculate the second invariant of deviatoric strain rate and stress       
    ε̇.II.c .= sqrt.((1/2) * (ε̇.xx.c.^2 + ε̇.yy.c.^2 + ε̇.zz.c.^2) + ε̇.xy.c.^2)
    τ.II.c .= sqrt.((1/2) * (τ.xx.c.^2 + τ.yy.c.^2 + τ.zz.c.^2) + τ.xy.c.^2)


    # Visualise solution fields
    p1 = heatmap(x.v, y.i, V.I.i',  aspect_ratio=1, title="V"  , xlims=extrema(x.v))
    p2 = heatmap(x.c, y.c, ε̇.II.c', aspect_ratio=1, title="ε̇II", xlims=extrema(x.v))
    p3 = heatmap(x.c, y.c, τ.II.c', aspect_ratio=1, title="τII", xlims=extrema(x.v))
    p4 = heatmap(x.c, y.c, P',      aspect_ratio=1, title="P"  , xlims=extrema(x.v))
    display(plot(p1, p2, p3, p4))
end

Mechanics2D()