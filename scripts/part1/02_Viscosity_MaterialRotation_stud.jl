using Plots, LinearAlgebra

ky = 1e3*365.25*24*3600

# @TODO: code the function that computes the second invariant
# assume a is a 3-component vector: [axx; ayy; axy]
# invII(a) = ...

function IsotropicViscousRheology(ε̇_cart, η, θ)

    # @TODO: Define the transformation matrix
    # T     = ...

    # Deviatoric strain rate in material coordinates 
    # @TODO: forward transformation
    # ε̇_mat    = ...

    # @TODO: Define the constitutive operator for isotropic viscous flow 
    # 𝐃_mat    = ...

    # @TODO: Evaluate deviatoric stress
    # τ_mat = ...

    # Deviatoric stress in Cartesian coordinates
    # @TODO: backward transformation
    # τ_cart   = ...

    return τ_cart, τ_mat
end

function (@main)()

        # Viscosity   
        η    = 5e22

        # Background strain rate
        ε̇bg  = 1e-15
    
        # Velocity gradient tensor
        L = [ε̇bg 0.; 0. -ε̇bg]

        # Deviatoric strain rate tensor in tensor form
        ϵ̇    = 1/2*(L + L') - 1/3*tr(L)*I(2)

        # Deviatoric strain rate tensor in Voigt form
        ε̇    = [ϵ̇[1,1]; ϵ̇[2,2]; ϵ̇[1,2]]

        # Define material orientations
        nθ         = 20
        θv         = LinRange(0, 360, nθ)
        τ_mat_all  = (
            xx  = zeros(nθ),
            xy  = zeros(nθ),
            II  = zeros(nθ),
        )

        # Loop on material orientations
        for iθ=1:nθ

            # Get angle
            θ  = θv[iθ] 

            # Compute stress 
            τ_cart, τ_mat = IsotropicViscousRheology(ε̇, η, θ)

            # Store data
            τ_mat_all.xx[iθ] = τ_mat[1]
            τ_mat_all.xy[iθ] = τ_mat[3]
            τ_mat_all.II[iθ] = invII(τ_mat)
        end

        # Plot stress state in material coordinates
        p1 = plot(xlabel="τxx_mat (MPa)", ylabel="τxy_mat (MPa)", aspect_ratio=1, xlims=(-120, 120), legend=:outerright)
        p1 = plot!(τ_mat_all.II .* cosd.(θv) ./ 1e6, τ_mat_all.II .* sind.(θv) ./ 1e6, label="Flow enveloppe" )
        p1 = scatter!(τ_mat_all.xx ./ 1e6, τ_mat_all.xy ./ 1e6, label="Stress state")

        # Plot invariant versus fabric orientation 
        p2 = plot(xlabel="θ (ᵒ)", ylabel="τii_mat (MPa)",  ylims=(0, 120))
        p2 = scatter!( θv, τ_mat_all.II./1e6, label="True invariant")

        plot(p1, p2, layout=(2,1))
end

main()

