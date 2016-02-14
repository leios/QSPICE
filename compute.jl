#=-------------compute.jl------------------------------------------------------#
#
# Purpose: To perform the computations necessary for each gate in the QSPICE
#          program
#
#-----------------------------------------------------------------------------=#

# First, the qubit type. Defined on Bloch sphere
type bloch
    r::Float64
    phi::Float64
    theta::Float64
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

# Implementing the not gate
function not(bit::qubit)
    N_gate = [Complex(0) Complex(1); Complex(1) Complex(0)]
    bit.element = N_gate * bit.element
    return bit
end

# Implementing a rotation gate
function rotation(bit::qubit, theta)
    R_gate = [Complex(1) Complex(0); Complex(0) Complex(exp(im * theta))]
    bit.element = R_gate * bit.element
    return bit
end

#=-----------------------------------------------------------------------------#
# MAIN
#-----------------------------------------------------------------------------=#

bit_bloch = bloch(1,0,pi)
bit = bloch_to_qubit(bit_bloch)
println(bit.element[1], '\t', bit.element[2])
bit = rotation(bit, pi)
println(bit.element[1], '\t', bit.element[2])

# List of functions that need to be implemented
#=
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

function measure()
end
=#
