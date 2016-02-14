#=-------------compute.jl------------------------------------------------------#
#
# Purpose: To perform the computations necessary for each gate in the QSPICE
#          program
#
#-----------------------------------------------------------------------------=#

# First, the qubit type. Defined on Bloch sphere
type bloch
    r::Float64
    theta::Float64
    phi::Float64
end

# Now to hold the information from the qubit
type qubit
    element::Array{Complex}
end

# Converting qubit type for manipulation later
# This should convert a qubit to the surface of the Bloch sphere
function bloch_to_qubit(angle::bloch)
    # Creating the cubit
    Bob = initialize()

    # define ray
    # Note: Not sure where to put these
    Bob.element[1] = cos(angle.theta / 2)
    Bob.element[2] = exp(im*angle.phi) * sin(angle.theta / 2)

    return Bob
end

function initialize()
    Bob = qubit([0,1])
    return Bob
end

Bob_Bloch = bloch(1,0.5 * pi,0.5 * pi)
Bob = bloch_to_qubit(Bob_Bloch)
println(Bob.element[1], '\t', Bob.element[2])

# List of functions that need to be implemented
#=
function hadamard()
end

function pauli_x()
end

function pauli_y()
end

function pauli_z()
end

function swap()
end

function xor()
end

function rotation()
end

function measure()
end
=#
