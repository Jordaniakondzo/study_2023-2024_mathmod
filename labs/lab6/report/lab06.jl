using Plots
using DifferentialEquations

a = 0.01
b = 0.02
N = 8439
y0 = 86
z0 = 25
x0 = N - y0 - z0

function ode_fn(du, u, p, t)
    x, y, z = u
    du[1] = -b*u[1]
    du[2] = a*u[1] - b*u[2]
    du[3] = b*u[1]
end

u0 = [x0, y0, z0]
tspan = (0.0, 200.0)
prob = ODEProblem(ode_fn, u0, tspan)
sol = solve(prob, dtmax=0.01)

X = [u[1] for u in sol.u]
Y = [u[2] for u in sol.u]
Z = [u[3] for u in sol.u]
T = [t for t in sol.t]

plt = 
    plot(
        layout=(1,2),
        dpi=300,
        legend=false)
    plot!(
        plt[1],
        T,
        X,
        title="решение уравнения S",
        color=:blue)
    plot!(
        plt[2],
        T,
        Y,
        label="решение уравнения I",
        color=:yellow)
    plot!(
        plt[2],
        T,
        Z,
        label="решение уравнения R",
        color=:red)

        savefig("lab6-1.png")