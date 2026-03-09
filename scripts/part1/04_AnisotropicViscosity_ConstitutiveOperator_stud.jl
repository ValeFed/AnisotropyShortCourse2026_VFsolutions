using Plots, LinearAlgebra, ForwardDiff, StaticArrays

ky = 1e3*365.25*24*3600

# @TODO: code the function that computes the second invariant
# assume a is a 3-component vector: [axx; ayy; axy]
# invII(a) = ...

# @TODO: define the anisotropic invariant
# invII_aniso(a, δ) = ...

function AnisotropicViscousRheology(ε̇_cart, η, θ, δ)

    # @TODO: Define the transformation matrix
    # T     = ...

    # Deviatoric strain rate in material coordinates 
    # @TODO: forward transformation
    # ε̇_mat    = ...

    # @TODO: Define the constitutive operator for transverse isotropic viscous flow 
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

        # Reduction of shear viscosity along weak direction
        δ    = 2

        # Background strain rate
        ε̇bg  = 1e-15
    
        # Velocity gradient tensor
        L = [ε̇bg 0.; 0. -ε̇bg]

        # Deviatoric strain rate tensor in tensor form
        ϵ̇    = 1/2*(L + L') - 1/3*tr(L)*I(2)

        # Deviatoric strain rate tensor in Voigt form
        ε̇    = [ϵ̇[1,1]; ϵ̇[2,2]; ϵ̇[1,2]]

        # Define material orientations
        nθ       = 50
        θv       = LinRange(0, 360, nθ)
        τii_cart = zeros(nθ)
        D11_mat  = zeros(nθ)
        D12_mat  = zeros(nθ)
        D13_mat  = zeros(nθ)
        D33_mat  = zeros(nθ)

        # Loop on material orientations
        for iθ=1:nθ

            # Get angle
            θ  = θv[iθ] 

            # Compute the stress in Cartesian coordinates
            τ_cart, τ_mat = AnisotropicViscousRheology(ε̇, η, θ, δ)

            # @TODO: Compute the rheological tensor in Cartesian coordinates using automatic differentiation
            # 𝐃_cart = ...

            # Display tensor
            display(𝐃_cart)

            # Store data
            D11_mat[iθ]  = 𝐃_cart[1,1]
            D12_mat[iθ]  = 𝐃_cart[1,2] 
            D13_mat[iθ]  = 𝐃_cart[1,3] 
            D33_mat[iθ]  = 𝐃_cart[3,3]
            τii_cart[iθ] = invII(τ_cart)

        end

        # Plot stress in material coordinates
        p1 = plot(xlabel="θ (ᵒ)", ylabel="D ij (Pa.s)")
        p1 = plot!(θv, D11_mat, label="D11")
        p1 = plot!(θv, D12_mat, label="D12")
        p1 = plot!(θv, D13_mat, label="D13")
        p1 = plot!(θv, D33_mat, label="D33")

        # Plot invariant versus fabric orientation 
        p2 = plot(xlabel="θ (ᵒ)", ylabel="Invariants (MPa)",  ylims=(0, 120))
        p2 = scatter!( θv, τii_cart ./1e6, label="True invariant")

        plot(p1, p2, layout=(2,1))
end

main()

