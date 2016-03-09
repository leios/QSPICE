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

using Pipe

import Base.copy
import Base.getindex
import Base.isapprox

# Macros
# Static variable macro, as described here:
# https://github.com/JuliaLang/julia/issues/15056
macro static(init)
  var = gensym()
  eval(current_module(), :(const $var = $init))
  var = esc(var)
  quote
    global $var
    $var
  end
end

# First, the Qubit type. Defined on Bloch sphere
type Bloch
    r::Float64
    theta::Float64
    phi::Float64
end

# Now to move the information from the Bloch to qubit
type Qubit
    vector::Vector{Complex{Float64}}
    bits::Int
end

# Initializes qubit with a particular state
function qubit(vector::Vector)
    bits = Int(log2(size(vector, 1)))
    return Qubit(vector, bits)
end

# Defining input for multiple states at once
qubit(s::Array...) = qubit(kron(s...))

# Constructs a qubit from its Bloch sphere representation
function qubit(angle::Bloch)
    vector::Vector{Complex} = ([cos(angle.theta / 2), 
                                exp(im*angle.phi) * sin(angle.theta / 2)])
    return qubit(vector)
end

function copy(q::Qubit)
    return Qubit(copy(q.vector), q.bits)
end

function getindex(qubit::Qubit, i::Int)
    return qubit.vector[i]
end

function isapprox(q1::Qubit, q2::Qubit, atol=0.0001)
    for (x, y) in zip(q1.vector, q2.vector)
        if !isapprox(x, y, atol=0.0001)
            return false
        end
    end
    return true
end

#=-----------------------------------------------------------------------------#
# FUNCTION
#-----------------------------------------------------------------------------=#
const IDENTITY = [1.0 + 0im 0; 0 1]

# Implement Hadamard gate
function hadamard(q::Qubit, bit::Int = 1)
    HADAMARD = @static (1 / sqrt(2)) * ([1.0 + 0im 1; 1 -1])

    # If there's exactly one bit, ignore the position and return 
    # the correct result
    if q.bits == 1
        return Qubit(HADAMARD * q.vector, q.bits)
    end

    gate_array::Array{Any, 1} = fill(IDENTITY, q.bits)
    gate_array[bit] = HADAMARD
    return Qubit(kron(gate_array...) * q.vector, q.bits)
end

# adding in rotations along the Bloch sphere

# Implementing the not gate
function not(q::Qubit)
    NOT = @static [0.0 + 0im 1; 1 0]
    return Qubit(NOT * q.vector, q.bits)
end

pauli_x(q::Qubit) = not(q)

function pauli_y(q::Qubit)
    P_Y = @static [0.0 0-1im; 0-1im 0]
    return Qubit(P_Y * q.vector, q.bits)
end

function pauli_z(q::Qubit)
    P_Z = @static [1.0 + 0im 0; 0 -1]
    return Qubit(P_Z * q.vector, q.bits)
end

# Implementing a rotation gate
function rotation(q::Qubit, theta::Real)
    rotation = [1.0 + 0im 0; 0 exp(im * theta)]
    return Qubit(rotation * q.vector, q.bits)
end

# Multi-qubit operations

# Implementing the swap gate
# Mutating version of swap, using the Julia convention of appending ! to 
# the mutating version
function swap!(q::Qubit, from::Int = 1, to::Int = 2)
    SWAP = @static [1.0 + 0im 0 0 0;
                    0 0 1 0;
                    0 1 0 0;
                    0 0 0 1]

    # If we want to swap on the same location we have nothing to do
    if from == to
        return
    end

    # In case there's only a single bit there's also nothing to swap
    if q.bits == 1
        return
    elseif q.bits == 2
        # If there are exactly 2 bits ignore the from and to 
        # arguments and swap them
        # alternatively we could throw an error here if from != 1 || to != 2
        q.vector = SWAP * q.vector
        return
    end

    # Normalize the swap interval, in case the user input is reversed
    low = min(from, to)
    high = max(from, to) - 1

    # The final swap gate that first swaps the lower bit into location then 
    # proceeds to restore the original order, while swapping the high bit down 
    # to its new location
    gate_array = fill(IDENTITY, q.bits - 1)
    gate_array[low] = SWAP
    swap_chain = kron(gate_array...)

    # Swap the lower bit up towards the high one
    for i = low + 1 : high
        gate_array = fill(IDENTITY, q.bits - 1)
        gate_array[i] = SWAP
        swap_chain = swap_chain * kron(gate_array...)
    end

    # Swap the high bit down to the original position of the low one
    for i = high - 1 : -1 : low
        gate_array = fill(IDENTITY, q.bits - 1)
        gate_array[i] = SWAP
        swap_chain = swap_chain * kron(gate_array...)
    end

    q.vector = swap_chain * q.vector
end

function swap(q::Qubit, from::Int = 1, to::Int = 2)
    ret = copy(q)
    swap!(ret, from, to)
    return ret
end

# Implementing the Cnot function
# Have not tested all cases, will look into later
# Only changed this to work with the rewrite, didn't improve the algorithm
function cnot(q::Qubit, control::Int, flip::Int)
    CNOT = @static [1.0 + 0im 0 0 0;
                    0 1 0 0;
                    0 0 0 1;
                    0 0 1 0]

    ret = copy(q)
    if ret.bits == 1
        return ret
    end

    cnot_chain = []

    # Swap the bits into the appropriate positions of 1 and 2 for
    # control and flip, respectively.

    # first swap everything into position
    swap!(ret, 1, control)
    if flip != 1
        swap!(ret, 2, flip)
    else
        swap!(ret, 2, control)
    end

    # performs cnot
    gate_array::Array{Any, 1} = fill(IDENTITY, q.bits - 1)
    gate_array[1] = CNOT

    ret.vector = kron(gate_array...) * ret.vector

    # swaps back
    if flip != 1
        swap!(ret, 2, flip)
    else
        swap!(ret, 2, control)
    end

    swap!(ret, 1, control)
    return ret
end

#=-----------------------------------------------------------------------------#
# MAIN
#-----------------------------------------------------------------------------=#

# Testing of single bit
println("Testing a single bit")
bit_bloch = Bloch(1,0,pi)
bit = qubit(bit_bloch)
@show bit
bit = rotation(bit, pi)
println("After rotation")
@show bit
println()

# Testing multi-Qubit operations
bit1 = [0;1]
bit2 = [1;0]
bit3 = [1;1]

# 2 qubit operations
println("Testing multi-qubit operations")
superbit = qubit(kron(bit1, bit2))
@show superbit

# We can chain operations using the @pipe macro. Use _ in place of the quantum bit.
# The operations are in order, from left to right
chaining = @pipe superbit |> hadamard(_, 1) |> swap(_, 1, 2) |> hadamard(_, 2)
@show chaining

println("\nTesting swap gate")
superbit3 = qubit(kron(bit1, bit2, bit3))
swapped3 = @pipe superbit3 |> swap(_, 1, 3)
@show superbit3
@show swapped3
println("Double swap equals original: ", isapprox(superbit3, swap(swapped3, 1, 3)))

println("\nTesting cnot gate")
cnot3 = cnot(superbit3, 2, 1)
@show superbit3
@show cnot3

# List of functions that need to be implemented
#=
function measure()
end
=#
