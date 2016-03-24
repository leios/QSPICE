#=-------------netlist_parser--------------------------------------------------#
#
# Purpose: To parse the netlist for the QSPICE program
#
# Provided by Gustorn:
# digit  ::= "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
# number ::= digit+
#
# parameter-list  ::= number ("," number)*
#
# gate-long  ::= "qubit" | "swap" | "not" | "cnot" | 
#                "rotate" | "hadamard" | "measure"
# gate-short ::= "q" | "s" | "n" | "cn" | "r" | "h" | "m"
# gate       ::= gate-long | gate-short
#
# gate-params      ::= "[" parameter-list "]"
# gate-connections ::= "(" parameter-list ")"
#
# element ::= number ":" gate gate-params? gate-connections?
#
# Human readable:
# <index>:<gate>[<optional parameters>](<connections>)
#
# 0:qubit[1,0,0,1]() 1:hadamard[2](0) 2:swap[1,2](0) 3:measure[](2)
#
#-----------------------------------------------------------------------------=#

include("compute.jl")

const FUNCTION_MAP = Dict{ASCIIString, Function}(
    "identity"       => qidentity,
    "hadamard"       => hadamard,
    "not"            => not,
    "cnot"           => cnot,
    "paulix"         => pauli_x,
    "pauliy"         => pauli_y,
    "pauliz"         => pauli_z,
    "swap"           => swap,
    "measure"        => measure,
    "qubit"          => qubit
)


# Creates a list of actions to execute from the netlist
function parse_nl(netlist_file::AbstractString)
    netlist = open(netlist_file, "r")
    commands = Dict()
    for line in readlines(netlist)
        terms = split(line)
        for term in terms

            println(term)
            # find action
            action = split(term, ":")

            # find arguments
            name, rest = split(action[2], "[")
            args, rest = split(rest, "]")

            if name != "qubit"
                args = map(x -> parse(Int, x), split(args,",",keep=false))
            else
                args = map(x -> parse(Float64, x), 
                           split(args,",",keep=false))
            end

            # find connections, -1 means it's not connected to anything
            conn = rest[2:end-1]
            if conn != ""
                conn = parse(Int, conn)
            else
                conn = -1
            end

            commands[parse(Int,action[1])] = Any[name, args, conn]
        end
    end
    println(commands)

    return commands
end

# Runs a single gate
function execute(commands)

# Maps all of the strings to functions

    temp_qubit = qubit([1,0,0,0])

    for i = 0:3
        println(i, '\t', temp_qubit)
        if commands[i][1] == "qubit"
            temp_qubit = FUNCTION_MAP[commands[i][1]](commands[i][2])
        elseif commands[i][1] == "measure"
            measurement = FUNCTION_MAP[commands[i][1]](temp_qubit,
                                                      commands[i][2]...)
            println(measurement)

        else
            temp_qubit = FUNCTION_MAP[commands[i][1]](temp_qubit,
                                                      commands[i][2]...)
        end
    end

end

# Implements list of commands
#execute_nl(commands)
#end

#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#

command_list = parse_nl("netlist.qnl")
execute(command_list)
