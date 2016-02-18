#=-------------compute.jl------------------------------------------------------#
#
# Purpose: To perform the computations necessary for each gate in the QSPICE
#          program
#
#   Notes: Multiple qubit states can be created with the kron function
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

# Define some gates (for testing)
type quantum_gate
    H
    N
    Id
    Swap
    Cont_n
end

H_gate = (1 / sqrt(2)) * [Complex(1) Complex(1); Complex(1) Complex(-1)]

Identity = [1 0; 0 1]

N_gate = [Complex(0) Complex(1); Complex(1) Complex(0)]

swap_gate = [Complex(1) Complex(0) Complex(0) Complex(0);
             Complex(0) Complex(0) Complex(1) Complex(0);
             Complex(0) Complex(1) Complex(0) Complex(0);
             Complex(0) Complex(0) Complex(0) Complex(1)]

Control_n = [Complex(1) Complex(0) Complex(0) Complex(0);
             Complex(0) Complex(1) Complex(0) Complex(0);
             Complex(0) Complex(0) Complex(0) Complex(1);
             Complex(0) Complex(0) Complex(1) Complex(0)]

gate = quantum_gate(H_gate, N_gate, Identity, swap_gate, Control_n)

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
    bit = qubit([1,0])
    return bit
end

# Implement Hadamard gate
function hadamard(bit::qubit, gate, num)
    # Create initial gate structure
    println(size(bit.element, 1))
    root_num = log2(size(bit.element, 1))
    gate_array = [Identity for i = 1:Int(root_num)]

    gate_array[num] = gate.H
 
    bit.element = kron(gate_array...) * bit.element

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

# Implementing the swap gate
function swap(bit::qubit)
    # Implementation of swap gate
    swap_gate = [Complex(1) Complex(0) Complex(0) Complex(0);
                 Complex(0) Complex(0) Complex(1) Complex(0);
                 Complex(0) Complex(1) Complex(0) Complex(0);
                 Complex(0) Complex(0) Complex(0) Complex(1)]

    bit.element = swap_gate * bit.element
end

#=-----------------------------------------------------------------------------#
# MAIN
#-----------------------------------------------------------------------------=#

# Testing of single bit

bit_bloch = bloch(1,0,pi)
bit = bloch_to_qubit(bit_bloch)
println(bit.element[1], '\t', bit.element[2])
bit = rotation(bit, pi)
println(bit.element[1], '\t', bit.element[2])

# Testing multi-qubit operations
bit1 = [0;1]
bit2 = [1;0]
bit3 = [1;0]

superbit = kron(bit1, bit2)
println("superbit is: ", superbit)

# H(1) -> swap -> H(2)
# We can chain operations together:
temp_bit = kron(Identity,H_gate) * swap_gate * kron(H_gate,Identity) * superbit

println(temp_bit)

# 3 qubit is similar to above:
superbit3 = qubit(kron(bit1, bit2, bit3))
superbit3 = hadamard(superbit3, gate, 2)
println(superbit3.element)

# List of functions that need to be implemented
#=
function pauli_x()
end

function pauli_y()
end

function pauli_z()
end

function measure()
end
=#
