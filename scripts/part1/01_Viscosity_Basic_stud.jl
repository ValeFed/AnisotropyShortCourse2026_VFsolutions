using Plots, LinearAlgebra

ky = 1e3*365.25*24*3600

# @TODO: code the function that computes the second invariant
# assume a is a 3-component vector: [axx; ayy; axy]
invII(a) = sqrt(0.5 * a[1]^2 + 0.5 * a[2]^2 + a[3]^2)

function (@main)()

    # Viscosity   
    η    = 1e23

    # Background strain rate
    ε̇bg  = 1e-15

    # Velocity gradient tensor
    L = [ε̇bg   0.; 
         0.   -ε̇bg] # Pure shear - only normal components, horizontal extension, vertical compression

    # @TODO: Define the deviatoric strain rate tensor
    ϵ̇   = 0.5 * (L + L') - (1/3) * tr(L) * I(2)


    # Deviatoric strain rate vector (Voigt notation)
    ε̇   = [ϵ̇[1,1]; ϵ̇[2,2]; ϵ̇[1,2]]

    # @TODO: Define the constitutive operator for isotropic viscous flow 
    D   = 2η*I(3)

    # @TODO: Evaluate deviatoric stress
    τ = D * ε̇

    # @TODO: Compute the deviatoric stress invariant
    τii = invII(τ)

    # Display deviatoric stress components
    display(τ)

    # Display deviatoric stress invariant
    display(τii)

    # Check
    display(τii == 2η * ε̇bg)
end

main()

