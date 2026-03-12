using Plots, LinearAlgebra

ky = 1e3*365.25*24*3600

# 1) code the function that computes the second invariant
# assume a is a 3-component vector: [axx; ayy; axy]
# invII(a) = ...

# invII_aniso(a, δ) = ...

function AnisotropicViscoelasticRheology(ε̇_cart, τ0_cart, η, G, θ, δ, Δt)

    # Transformation matrix
    c, s  = cosd(θ), sind(θ)
    # T     = ...

    # @TODO: Deviatoric strain rate in material coordinates (forward transformation)
    # ε̇_mat     = ...

    # @TODO: Old deviatoric stress in material coordinates (forward transformation)
    # τ0_mat    = ...

    # @TODO: Effective strain rate accounting for old deviatoric stress
    # ε̇_mat_eff = ...

    # @TODO: Effective viscosity
    # η_eff = ...

    # @TODO: Constitutive matrix in material coordinates
    # 𝐃_mat    = ...

    # @TODO: Deviatoric stress in material coordinates
    # τ_mat = ...

    # @TODO: Deviatoric stress in Cartesian coordinates (backward transformation)
    # τ_cart   =...

    return τ_cart, τ_mat
end


function (@main)()


        # Number of time steps
        nt   = 100

        # Time step
        Δt   = 1e11

        # Viscosity   
        η    = 5e22
        
        # Shear modulus
        G    = 3e10

        # Reduction of shear viscosity/modulus along weak direction
        δ    = 2

        # Background strain rate
        ε̇bg  = 1e-15
    
        # Velocity gradient tensor
        L = [ε̇bg 0.; 0. -ε̇bg]

        # Deviatoric strain rate tensor in tensor form
        ϵ̇    = 1/2*(L + L') - 1/3*tr(L)*I(2)

        # Deviatoric strain rate tensor in Voigt form
        ε̇    = [ϵ̇[1,1]; ϵ̇[2,2]; ϵ̇[1,2]]

        # Initial deviatoric stress
        τ_cart  = zeros(3)
        τ0_cart = zeros(3)

        # Define material orientations
        nθ       = 50
        θv       = LinRange(0, 360, nθ)
        τxx_mat  = zeros(nθ, nt)
        τxy_mat  = zeros(nθ, nt)
        τii_cart = zeros(nθ, nt)
        τii_mat  = zeros(nθ, nt)

        # Loop on material orientations
        for iθ=1:nθ

            # Get angle
            θ  = θv[iθ] 

            # Set initial condition
            τ_cart .= 0.0

            # Time loop
            for it=1:nt

                # From previous time step
                τ0_cart .= τ_cart

                # Compute the stress in Cartesian coordinates
                τ_cart, τ_mat = AnisotropicViscoelasticRheology(ε̇, τ0_cart, η, G, θ, δ, Δt)
                
                # Store data
                τxx_mat[iθ, it]  = τ_mat[1]
                τxy_mat[iθ, it]  = τ_mat[3]
                τii_cart[iθ, it] = invII(τ_cart)
                τii_mat[iθ, it]  = invII_aniso(τ_mat, δ)
            end
            
        end

        # Plot stress in material coordinates
        p1 = plot(xlabel="τxx_mat (MPa)", ylabel="τxy_mat (MPa)", aspect_ratio=1, xlims=(-120, 120), legend=:outerright)
        p1 = plot!(τii_mat  .* cosd.(θv) ./ 1e6, τii_mat/δ .* sind.(θv) ./ 1e6, label=:none)#, label="Flow enveloppe (material invariant)")
        p1 = scatter!(τxx_mat ./ 1e6, τxy_mat ./ 1e6, markersize=2, label=:none)#, label="Stress state")

        # Plot invariant versus fabric orientation at final time
        p2 = plot(xlabel="θ (ᵒ)", ylabel="τII (MPa)",  ylims=(0, 120))
        p2 = scatter!( θv, τii_cart[:,end] ./1e6, label="True invariant")
        p2 = scatter!( θv, τii_mat[:,end]  ./1e6, label="Material pseudo-invariant")

        # Plot invariant versus time 
        ind45 = argmin( abs.(θv .- 45))
        p3 = plot(xlabel="t (ky)", ylabel="τII (MPa)",  ylims=(0, 120))
        for i = 1:ind45
            p3 = plot!((1:nt)*Δt / ky, τii_cart[i, :] ./1e6, label=:none, linestyle=:dash )
        end
        p3 = plot!((1:nt)*Δt / ky, τii_cart[1,     :] ./1e6, label="Strong direction", linewidth=2 )
        p3 = plot!((1:nt)*Δt / ky, τii_cart[ind45, :] ./1e6, label="Weak direction", linewidth=2 )
        plot(p1, p2, p3, layout=(3,1))
end

main()

