using PauliPropagation
using Plots

Alist = []
# nqubits = 5
nlayers = 15
range = 4:10
alpha = 1.5
for nqubits in range
    println("nqubits = ", nqubits)
    ## define the observable
    # here Z......I
    local observable = PauliSum(nqubits)
    add!(observable, [:Y,:Z,:X,:X], [1,2,3,4], 1.0)
    # add!(observable, [:X, :Y, :Z], [1, 2, 3], 1.3)
    # local nlayers = 50
    local topology = bricklayertopology(nqubits; periodic=true)
    local circuit =su4circuit(nqubits, nlayers; topology=topology)
    local dt = 3.14159265358979323846
    local nparams = countparameters(circuit)
    local parameters = rand(nparams) * dt
    local pauli_sum = propagate(circuit, observable, parameters)
    # 提取 Float64 部分
    println(pauli_sum)
    local my_values = [x^2 for x in collect(values(pauli_sum.terms))]
    local my_length = length(my_values)
    local sorted_my_values = sort(my_values, rev=true)
    local A_values = [(i^alpha)*sorted_my_values[i] for i in 1:my_length]
    # 从大到小排序
    local sorted_A_values = sort(A_values, rev=true)
    # 绘制图形
    # p = plot(sorted_my_values, seriestype=:scatter, title="Sorted Pauli Terms", xlabel="Index", ylabel="Value")
    # p2 = plot(A_values, seriestype=:scatter, title="A Terms", xlabel="Index", ylabel="Value")
    local A = sorted_A_values[1]
    println("A = ", A)
    inverse_function(x, A_cons) = A_cons ./ x^alpha
    xdata = 1:my_length
    # plot!(xdata, inverse_function(xdata, A), label="A / x")
    # display(p)
    # display(p2)
    println("my_length=", my_length)
    println("sum_of_values=", sum(my_values))
    push!(Alist, A)
end
# 绘制 A vs nqubits 的图形
scatter(range, Alist, title="A vs nqubits", xlabel="nqubits", ylabel="A", label="Data")
# plot(range, Alist, title="A vs nlayers", xlabel="nlayers", ylabel="A", label="Data")
