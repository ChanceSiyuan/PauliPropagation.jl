# using Pkg
# Pkg.add("Bits")
# Pkg.add("Plots")
using PauliPropagation
using Plots
using Random

Random.seed!(1234)
nqubits = 5
alpha = 1.5
observable = PauliSum(nqubits)
add!(observable, [:Y,:Z,:X], [1,2,3], 1.0)
nlayers = 50
topology = bricklayertopology(nqubits; periodic=true)
circuit =su4circuit(nqubits, nlayers; topology=topology)
dt = 3.14159265358979323846
nparams = countparameters(circuit)
parameters = rand(nparams) * dt
pauli_sum = propagate(circuit, observable, parameters)
println(pauli_sum)
my_values = [x^2 for x in collect(values(pauli_sum.terms))]
println("sum_of_values=", sum(my_values))
my_length = length(my_values)
println("my_length=", my_length)
sorted_my_values = sort(my_values, rev=true)
# range = 1:1
# Alist = []
# for r in range
#     A_values = [(i^r)*sorted_my_values[i] for i in 1:my_length]
#     println(A_values)
#     # 从大到小排序
#     sorted_A_values = sort(A_values, rev=true)
#     println(sorted_A_values)
#     A = sorted_A_values[1]
#     push!(Alist, A)
# end
# # 绘制 A vs r 的图形
# scatter(range, Alist, title="A vs r", xlabel="r", ylabel="A", label="Data",yscale=:log10)
# # plot(range, Alist, title="A vs nlayers", xlabel="nlayers", ylabel="A", label="Data")
# Plot for r = 0
r = 0
A_values = [(i^r)*sorted_my_values[i] for i in 1:my_length]
global pic = plot(1:my_length, A_values, title="function vs i", xlabel="i", ylabel="i^r * f", label="r = 0")

for r in 0.1:0.1:1.5
    local f = [(i^r)*sorted_my_values[i] for i in 1:my_length]
    plot!(pic,1:my_length, f, label="r =$(r)")
end
display(pic)
# # Plot for r = 2
# for r in 0.1:0.1
#     f = [(i^r)*sorted_my_values[i] for i in 1:my_length]
#     plot!(1:my_length, f, label="i^$(r) * f")
# end
