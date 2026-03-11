using Plots, LinearAlgebra

ky = 1e3*365.25*24*3600

# @TODO: code the function that computes the second invariant
# assume a is a 3-component vector: [axx; ayy; axy]
# invII(a) = ...

function (@main)()

    # Viscosity   
    η    = 1e23

    # Background strain rate
    ε̇bg  = 1e-15

    # Velocity gradient tensor
    L = [ε̇bg   0.; 
         0.   -ε̇bg]

    # @TODO: Define the deviatoric strain rate tensor
    # ϵ̇   = ...

    # Deviatoric strain rate vector (Voigt notation)
    ε̇   = [ϵ̇[1,1]; ϵ̇[2,2]; ϵ̇[1,2]]

    # @TODO: Define the constitutive operator for isotropic viscous flow 
    # 𝐃   = ...

    # @TODO: Evaluate deviatoric stress
    # τ = ...

    # @TODO: Compute the deviatoric stress invariant
    # τii = ...

    # Display deviatoric stress components
    display(τ)

    # Display deviatoric stress invariant
    display(τii)
end

main()

