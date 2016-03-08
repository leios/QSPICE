#=-------------test.jl---------------------------------------------------------#
#
# Purpose: To test the functions in the compute.jl file with the qubit
#          operations
#
#   Notes: FactCheck implements an assert-like function in Julia
#          Iterators allows for "for i in product(0:1, 0:1, 0:1)"
#
#-----------------------------------------------------------------------------=#

# Including works for in-directory files
include("compute.jl")

# Using works for libraries
using FactCheck
using Iterators

const QUBITS = Array[[0,0], [0,1], [1,0], [1,1]]

# Generates random qubit for later
function random_qubit(num_bits)
    return qubit(reduce(kron, [QUBITS[rand(1:length(QUBITS))] for _ in 1:num_bits]))
end

function random_qubit_list(num_bits, n)
    return [random_qubit(num_bits) for _ in 1:n]
end

# Generates all distinct quantum bit vectors for num_bits
function exhaustive_list(num_bits)
    arrays = ([(Array[x...]) for x in product([QUBITS for i=1:num_bits]...)])
    states = imap(x -> kron(x...), arrays) |> distinct |> x -> imap(qubit, x) |> collect
    return states
end

roughly(q1::Qubit) = (q2::Qubit) -> isapprox(q1, q2)

# Fact checking / asserting that the swaps happen as expected
facts("Swap gate tests") do
    context("Test for a single qubit") do
        q = qubit([0,1])
        @fact q --> roughly(swap(q))
        @fact q --> roughly(swap(swap(q)))
        @fact q --> roughly(swap(q, 4, 5))
        @fact q --> roughly(swap(swap(q, 1, 4), 4, 5))
    end

    context("Exhaustive test for 2 qubits") do
        qubits = exhaustive_list(2)
        for q in qubits
            @fact q --> roughly(swap(swap(q)))
        end
    end

    context("Exhaustive test for 3 qubits") do
        qubits = exhaustive_list(3)
        for q in qubits
            for j = 1:3, k = 1:3
                @fact q --> roughly(swap(swap(q, j, k), j, k))
            end
        end
    end

    context("Exhaustive test for 4 qubits") do
        qubits = exhaustive_list(4)
        for q in qubits
            for j = 1:4, k = 1:4
                @fact q --> roughly(swap(swap(q, j, k), j, k))
            end
        end
    end

    context("Exhaustive test for 5 qubits") do
        qubits = exhaustive_list(5)
        swap(swap(qubits[1], 1, 4), 1, 5)
        for q in qubits
            for j = 1:5, k = 1:5
                @fact q --> roughly(swap(swap(q, j, k), j, k))
            end
        end
    end
end
