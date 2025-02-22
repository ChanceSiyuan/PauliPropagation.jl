# using Pkg
# Pkg.add("Bits")
# Pkg.add("Plots")
using PauliPropagation
using Plots
using Random

function Plot_transition_point!(pic, circuit, observable, parameters; label="Data")
    #开始Pauli反向传播
    pauli_sum = propagate(circuit, observable, parameters)
    #提取最终Pauli项的系数值
    my_values = [x^2 for x in collect(values(pauli_sum.terms))]
    my_length = length(my_values)
    # println("sum_of_values=", sum(my_values)) #校验是否归一
    # println("my_length=", my_length)
    sorted_my_values = sort(my_values, rev=true)
    range = 0:0.01:0.5
    Alist = []
    for r in range
        local A_values = [(i^r)*sorted_my_values[i] for i in 1:my_length]
        sorted_A_values = sort(A_values, rev=true) # 从大到小排序
        A = sorted_A_values[1] #取最大值
        push!(Alist, A)
    end
    # 绘制 A vs r 的图形
    scatter!(pic, range, Alist, label=label)
end


nqubits = 10
##############################################################################
#生成随机一排2-qubit unitary拼成的线路
Random.seed!(1234)
nlayers = 1
topology = bricklayertopology(nqubits; periodic=true)
circuit1 = su4circuit(nqubits, nlayers; topology=topology)
dt = 3.14159265358979323846
nparams = countparameters(circuit1)
parameters1 = rand(nparams) * dt
#确定最终要测量的Pauli算符
observable1 = PauliSum(nqubits)
add!(observable1, [:X for i in 1:nqubits], [i for i in 1:nqubits], 1.0)
##############################################################################
#生成另一排2-qubit unitary拼成的线路
Random.seed!(5678)
circuit2 = su4circuit(nqubits, nlayers; topology=topology)
parameters2 = rand(nparams) * dt
#确定最终要测量的Pauli算符
observable2 = PauliSum(nqubits)
add!(observable2, [:X for i in 1:nqubits], [i for i in 1:nqubits], 1.0)
##############################################################################
#生成一排T gate的线路
circuit3::Vector{Gate} = []
append!(circuit3,(PauliRotation(:Z, i) for i in 1:nqubits))
parameters3 = [pi/4 for i in 1:nqubits]
#确定最终要测量的Pauli算符
observable3 = PauliSum(nqubits)
add!(observable3, [:X for i in 1:nqubits], [i for i in 1: nqubits], 1.0)
##############################################################################
# Create an empty plot
pic = plot(title="A vs r", xlabel="r", ylabel="A", yscale=:log10)
# Plot the first dataset
Plot_transition_point!(pic, circuit1, observable1, parameters1; label="random circuit 1")
# Plot the second dataset
Plot_transition_point!(pic, circuit2, observable2, parameters2; label="random circuit 2")
# Plot the Third dataset
Plot_transition_point!(pic, circuit3, observable3, parameters3; label="T-layer circuit")

# Display the plot
display(pic)



