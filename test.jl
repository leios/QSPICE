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

# Generates random qubit for later
function random_qubit(length)
    # initialize array to a random set of 1/0's
    r_array = [rand(0:1) for i in 1:length]

    # create / return qubit
    r_qubit = qubit(r_array)
    return r_qubit
end

# Fact checking / asserting that the swaps happen as expected
facts("Checking the swap gate for all cases up to qubit length 5.") do

    # Creating initial states
    state_1 = [0;1]
    state_0 = [1;0]
    
    # using iterator
    for i in product(0:1, 0:1, 0:1, 0:1, 0:1)

        # generating qubit states from iterator
        q_states = [state_0 for j in 1:length(i)]
        for j in i
            if j != 0
                q_states[j] = state_1
            end
        end

        test_qubit = qubit(kron(q_states...))

        println(i)
        for j = 1:5, k = 1:5
            println(j, '\t', k)
            @fact test_qubit --> swap(swap(test_qubit,gate,j,k),gate,j,k)
        end
    end

end
