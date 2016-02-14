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

#=-----------------------------------------------------------------------------#
# FUNCTION
#-----------------------------------------------------------------------------=#

# Converting qubit type for manipulation later
# This should convert a qubit to the surface of the Bloch sphere
function bloch_to_qubit(angle::bloch)
    # Creating the cubit
    bit = initialize()

    # define ray
    # Note: Not sure where to put these
    bit.element[1] = cos(angle.theta / 2)
    bit.element[2] = exp(im*angle.phi) * sin(angle.theta / 2)

    return bit
end

# initializes a single qubit in the |0> state
function initialize()
    bit = qubit([0,1])
    return bit
end

# Implement Hadamard gate
function hadamard(bit::qubit)
    # This is a predefined Hadamard gate
    H_gate = (1 / sqrt(2)) * [Complex(1) Complex(1); Complex(1) Complex(-1)]
    bit.element = H_gate * bit.element
    return bit
end

#=-----------------------------------------------------------------------------#
# MAIN
#-----------------------------------------------------------------------------=#

bit_bloch = bloch(1,0,0)
bit = bloch_to_qubit(bit_bloch)
bit = hadamard(bit)
println(bit.element[1], '\t', bit.element[2])

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
