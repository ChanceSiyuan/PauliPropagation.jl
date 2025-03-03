# using Pkg
# Pkg.add("Bits")
# Pkg.add("Plots")
# Pkg.add("LaTeXStrings")
using PauliPropagation
using Plots
using Random
using LaTeXStrings

function Plot_transition_point!(pic, circuit, observable, parameters; label="Data")
    #开始Pauli反向传播
    pauli_sum = propagate(circuit, observable, parameters)
    print(pauli_sum)
    #提取最终Pauli项的系数值
    my_values = [x^2 for x in collect(values(pauli_sum.terms))]
    my_length = length(my_values)
    # println("sum_of_values=", sum(my_values)) #校验是否归一
    # println("my_length=", my_length)
    sorted_my_values = sort(my_values, rev=true)
    range = -1:0.01:0.2
    Alist = []
    for r in range
        local A_values = [(i^(1 + r))*sorted_my_values[i] for i in 1:my_length]
        sorted_A_values = sort(A_values, rev=true) # 从大到小排序
        A = sorted_A_values[1] #取最大值
        push!(Alist, A)
    end
    # 绘制 A vs r 的图形
    plot!(pic, range, Alist, label=label)
end

function noisy_su4circuit(nqubits::Integer, nlayers::Integer; topology=nothing,noise_rate::Real,noise_type::Integer)
    
"""
noisy_su4circuit(nqubits::Integer, nlayers::Integer; topology=nothing)
Create a circuit that consists of layers of SU(4) gates on a given topology and depolarizing noise layers on each qubit. 
If no topology is specified, a bricklayer topology is used.

"""
    circuit::Vector{Gate} = []
    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end
    if noise_type == 1
        append!(circuit, (FrozenGate(DepolarizingNoise(qind), noise_rate) for qind in 1:nqubits))
        for nl in 1:nlayers
            for pair in topology
                appendSU4!(circuit, pair)
            end
            append!(circuit, (FrozenGate(DepolarizingNoise(qind), noise_rate)  for qind in 1:nqubits))
        end
    elseif noise_type == 2
        append!(circuit, (FrozenGate(DephasingNoise(qind), noise_rate) for qind in 1:nqubits))
        for nl in 1:nlayers
            for pair in topology
                appendSU4!(circuit, pair)
            end
        end
            append!(circuit, (FrozenGate(DephasingNoise(qind), noise_rate)  for qind in 1:nqubits))
    elseif noise_type == 3
        append!(circuit, (FrozenGate(AmplitudeDampingNoise(qind), noise_rate) for qind in 1:nqubits))
        for nl in 1:nlayers
            for pair in topology
                appendSU4!(circuit, pair)
            end
            append!(circuit, (FrozenGate(AmplitudeDampingNoise(qind), noise_rate)  for qind in 1:nqubits))
        end
    end

    return circuit
end

function noisy_Clifford_T_circuit(nqubits::Integer, nlayers::Integer; topology=nothing,noise_rate::Real)
    circuit::Vector{Gate} = []
    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end
    append!(circuit, (FrozenGate(DepolarizingNoise(qind), noise_rate)  for qind in 1:nqubits))
    for nl in 1:nlayers
        for pair in topology
            if rand(Bool)
                push!(circuit, CliffordGate(:CNOT, pair))
            end
        end
        for indx in 1:nqubits
            if rand(Bool)
                push!(circuit, CliffordGate(:H, indx))
            end
            if rand(Bool)
                push!(circuit, CliffordGate(:S, indx))
            end
        end
        append!(circuit, (FrozenGate(PauliRotation(:Z, i),pi/4) for i in 1:nqubits))
        append!(circuit, (FrozenGate(DepolarizingNoise(qind), noise_rate)  for qind in 1:nqubits))
    end
    return circuit
end

function one_layer_of_T_gate(nqubits::Integer;noise_rate::Real,noise_type::Integer)
    circuit::Vector{Gate} = []
    append!(circuit,(FrozenGate(PauliRotation(:Z, i),pi/4) for i in 1:nqubits))
    if noise_type == 1
        append!(circuit, (FrozenGate(DepolarizingNoise(qind), noise_rate)  for qind in 1:nqubits))
    elseif noise_type == 2
        append!(circuit, (FrozenGate(DephasingNoise(qind), noise_rate)  for qind in 1:nqubits))
    elseif noise_type == 3
        append!(circuit, (FrozenGate(AmplitudeDampingNoise(qind), noise_rate)  for qind in 1:nqubits))
    end
    return circuit
end



Random.seed!(1234)
global nqubits = 8
global topology = bricklayertopology(nqubits; periodic=true)
global dt = 3.14159265358979323846
global nlayers = 10
##############################################################################
# # Create an empty plot
# pic1 = plot(title="T+Depolarization noise", xlabel=L"\alpha", ylabel=L"C(\alpha)", yscale=:log10,legend=:outerright)
# #确定最终要测量的Pauli算符
# observable3 = PauliSum(nqubits)
# parameters3 =[]
# add!(observable3, [:X for i in 1:nqubits], [i for i in 1: nqubits], 1.0)
# for noise_rate in 0:0.02:0.2
#     #生成一排T gate的线路
#     circuit3=one_layer_of_T_gate(nqubits;noise_rate,noise_type = 1)
#     # # Plot the Third dataset
#     Plot_transition_point!(pic1, circuit3, observable3, parameters3; label="p= $noise_rate")
# end
# display(pic1)
##############################################################################
# # Create an empty plot
# pic2 = plot(title="T+Dephasing noise", xlabel=L"\alpha", ylabel=L"C(\alpha)", yscale=:log10,legend=:outerright)
# #确定最终要测量的Pauli算符
# observable3 = PauliSum(nqubits)
# parameters3 =[]
# add!(observable3, [:X for i in 1:nqubits], [i for i in 1: nqubits], 1.0)
# for noise_rate in 0:0.02:0.2
#     #生成一排T gate的线路
#     circuit3=one_layer_of_T_gate(nqubits;noise_rate,noise_type = 2)
#     # # Plot the Third dataset
#     Plot_transition_point!(pic2, circuit3, observable3, parameters3; label="p= $noise_rate")
# end
# display(pic2)
##############################################################################
# # Create an empty plot
# pic3 = plot(title="T+Depolarizing noise", xlabel=L"\alpha", ylabel=L"C(\alpha)", yscale=:log10,legend=:outerright)
# #确定最终要测量的Pauli算符
# observable3 = PauliSum(nqubits)
# parameters3 =[]
# add!(observable3, [:X for i in 1:nqubits], [i for i in 1: nqubits], 1.0)
# for noise_rate in 0:0.02:0.2
#     #生成一排T gate的线路
#     circuit3=one_layer_of_T_gate(nqubits;noise_rate,noise_type = 3)
#     # # Plot the Third dataset
#     Plot_transition_point!(pic3, circuit3, observable3, parameters3; label="p= $noise_rate")
# end
# display(pic3)
##############################################################################

# #生成一排Clifford+ T gate的线路
# circuit_cliffordT=noisy_Clifford_T_circuit(nqubits,nlayers;topology = nothing,noise_rate = 0.2)
# #确定最终要测量的Pauli算符
# observable_cliffordT = PauliSum(nqubits)
# add!(observable_cliffordT, [:X for i in 1:nqubits], [i for i in 1: nqubits], 1.0)
# parameters_cliffordT =[]
# # # Plot the Third dataset
# Plot_transition_point!(pic, circuit_cliffordT, observable_cliffordT, parameters_cliffordT; label="Clifford+T circuit")

##############################################################################
# #生成随机一排2-qubit unitary拼成的线路circuit1
# circuit1 = su4circuit(nqubits, nlayers; topology=topology)
# nparams = countparameters(circuit1)
# print(nparams)
# parameters1 = rand_para * dt
# #确定最终要测量的Pauli算符
# observable1 = PauliSum(nqubits)
# # add!(observable1, [:X for i in 1:nqubits], [i for i in 1:nqubits], 1.0)
# add!(observable1, [:X], [7], 1.0)
# # # Plot the first dataset
# Plot_transition_point!(pic, circuit1, observable1, parameters1; label="Noiseless circuit")
# ##############################################################################
#生成另一排含噪声2-qubit unitary拼成的线路circuit2
#改变noise的大小
pic = plot(title="su(4)+Depolarizing noise", xlabel=L"\alpha", ylabel=L"C(\alpha)", yscale=:log10,legend=:outerright)
for noise_rate in 0:0.01:0.1
    local circuit2 = noisy_su4circuit(nqubits, nlayers; topology=topology,noise_rate = noise_rate,noise_type = 3)
    local nparams = countparameters(circuit2)
    local rand_para = rand(nparams)
    local parameters2 = rand_para * dt
    #确定最终要测量的Pauli算符
    local observable2 = PauliSum(nqubits)
    # add!(observable2, [:X for i in 1:nqubits], [i for i in 1:nqubits], 1.0)
    add!(observable2, [:X], [7], 1.0)
    # # Plot the second dataset
    Plot_transition_point!(pic, circuit2, observable2, parameters2; label="Noisy circuit, p= $noise_rate")
end
display(pic)
# ############################################################################


# ##############################################################################
# #生成另一排含噪声2-qubit unitary拼成的线路circuit2
# #改变layer的大小
# for nlayers in 1:4:28
#     local circuit2 = noisy_su4circuit(nqubits, nlayers; topology=topology,noise_rate = 0.02)
#     local nparams = countparameters(circuit2)
#     local rand_para = rand(nparams)
#     local parameters2 = rand_para * dt
#     #确定最终要测量的Pauli算符
#     local observable2 = PauliSum(nqubits)
#     # add!(observable2, [:X for i in 1:nqubits], [i for i in 1:nqubits], 1.0)
#     add!(observable2, [:X], [7], 1.0)
#     # # Plot the second dataset
#     Plot_transition_point!(pic, circuit2, observable2, parameters2; label="Noisy circuit, layer number= $nlayers")
# end
# ############################################################################
# # Display the plot




