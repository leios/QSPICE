#=-------------compute.jl------------------------------------------------------#
#
# Purpose: To perform the computations necessary for each gate in the QSPICE
#          program
#
#   Notes: Multiple Qubit states can be created with the kron function
#          We can only swap between adjacent positions
#          This should just be swapping the xth and yth row in the ID matrix
#
#-----------------------------------------------------------------------------=#

import Base.getindex

# First, the Qubit type. Defined on Bloch sphere
type Bloch
    r::Float64
    phi::Float64
    theta::Float64
end

# Now to hold the information from the Qubit
type Qubit
    element::Array{Complex}
end

function qubit()
    return Qubit([0, 1])
end

function qubit(angle::Bloch)
    return Qubit([cos(angle.theta / 2),
                  exp(im*angle.phi) * sin(angle.theta / 2)])
end

function approx_eq(q1::Qubit, q2::Qubit)
    for (x, y) in zip(q1.element, q2.element)
        if abs2(x - y) > 0.0001
            return false
        end
    end
    return true
end

function getindex(qubit::Qubit, i::Int)
    return qubit.element[i]
end

const IDENTITY = [Complex(1) 0; 0 1]

#=-----------------------------------------------------------------------------#
# FUNCTION
#-----------------------------------------------------------------------------=#

# Implement Hadamard gate
function hadamard(bit::Qubit, num = 1)
    const hadamard_matrix = (1 / sqrt(2)) * [Complex(1) Complex(1); Complex(1) Complex(-1)]

    # Create initial gate structure
    root_num = Int(log2(size(bit.element, 1)))
    if root_num == 0
        return Qubit(hadamard_matrix * bit.element)
    end

    gate_array::Array{Any,1} = [IDENTITY for i = 1:root_num]
    gate_array[num] = hadamard_matrix
    return Qubit(kron(gate_array...) * bit.element)
end

# Implementing the not gate
function not(bit::Qubit)
    const not_matrix = [Complex(0) Complex(1); Complex(1) Complex(0)]
    return Qubit(not_matrix * bit.element)
end

# Implementing a rotation gate
function rotation(bit::Qubit, theta)
    const rotation_matrix = [Complex(1) Complex(0); Complex(0) Complex(exp(im * theta))]
    return Qubit(rotation_matrix * bit.element)
end

# Implementing the swap gate
function swap(bit::Qubit, from = 1, to = 2)
    const swap_matrix = [Complex(1) Complex(0) Complex(0) Complex(0);
                         Complex(0) Complex(0) Complex(1) Complex(0);
                         Complex(0) Complex(1) Complex(0) Complex(0);
                         Complex(0) Complex(0) Complex(0) Complex(1)]

    if from == to
        return Qubit(copy(bit.element))
    end

    # Create initial gate structure
    # Number of bits = root_num
    root_num = Int(log2(size(bit.element, 1)))

    if root_num <= 1
        return Qubit(copy(bit.element))
    elseif root_num == 2
        test = swap_matrix * bit.element
        return Qubit(swap_matrix * bit.element)
    end

    # moving from low to high
    low = min(from, to)
    high = max(from, to) - 1

    # final gate for multiplication
    swap_chain = []

    for i = low : high
        gate_array = [IDENTITY for i = 1 : root_num - 1]
        gate_array[i] = swap_matrix
        if i == low
            swap_chain = kron(gate_array...)
        else
            swap_chain = swap_chain * kron(gate_array...)
        end
    end

    for i = high - 1 : -1 : low
        gate_array = [gate.Id for i = 1:root_num - 1]
        gate_array[i] = gate.Swap
        swap_chain = swap_chain * kron(gate_array...)
    end

    return Qubit(swap_chain * bit.element)
end

# Implementing the Cnot function
# Have not tested all cases, will look into later
function cnot(bit::Qubit, control, flip)
    const cnot_matrix = [Complex(1) Complex(0) Complex(0) Complex(0);
                         Complex(0) Complex(1) Complex(0) Complex(0);
                         Complex(0) Complex(0) Complex(0) Complex(1);
                         Complex(0) Complex(0) Complex(1) Complex(0)]

    root_num = Int(log2(size(bit.element, 1)))

    if root_num <= 1
        return bit
    end

    cnot_chain = []

    # Swap the bits into the appropriate positions of 1 and 2 for
    # control and flip, respectively.

    # first swap everything into position
    bit = swap(bit, 1, control)
    if flip != 1
        bit = swap(bit, 2, flip)
    else
        bit = swap(bit, 2, control)
    end

    # performs cnot
    gate_array = [IDENTITY for i = 1:root_num-1]
    gate_array[1] = cnot_matrix

    bit.element = kron(gate_array...) * bit.element

    # swaps back
    if flip != 1
        bit = swap(bit, 2, flip)
    else
        bit = swap(bit, 2, control)
    end

    bit = swap(bit, gate, 1, control)
    return bit
end

#=-----------------------------------------------------------------------------#
# MAIN
#-----------------------------------------------------------------------------=#

# Testing of single bit

bit_bloch = Bloch(1,0,pi)
bit = qubit(bit_bloch)
println(bit[1], '\t', bit[2])
bit = rotation(bit, pi)
println(bit[1], '\t', bit[2])

# Testing multi-Qubit operations
bit1 = [0;1]
bit2 = [1;0]
bit3 = [0;1]

superbit = Qubit(kron(bit1, bit2))
println("superbit is: ", superbit)

# H(1) -> swap -> H(2)
# We can chain operations together (in reverse order):
# temp_bit = superbit |> q -> hadamard(q, 1) |> q -> swap(q, 1, 2) |> q -> hadamard(q, 2)

temp_bit = superbit |> swap
@show superbit
@show temp_bit

println(approx_eq(superbit, temp_bit))

# 3 Qubit is similar to above:
#superbit3 = Qubit(kron(bit1, bit2, bit3))
#=
println("original superbit")
println(superbit3.element)
superbit3 = swap(superbit3, gate, 1, 3)
println("superbit after swap")
println(superbit3.element)
=#

#println(superbit3 == swap(swap(superbit3, gate, 1, 3), gate, 1, 3))
#
#superbit3 = cnot(superbit3, gate, 2, 1)
#println("superbit after cnot")
#println(superbit3.element)

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
