using DifferentialEquations, Plots

# Определение функции модели хищник-жертва
function pred_prey!(du, u, p, t)
    du[1] = -0.13*u[1] + 0.04*u[1]*u[2] # изменение популяции жертв dx/dt
    du[2] = 0.31*u[2] - 0.042*u[1]*u[2] # изменение популяции хищников dy/dt
end

# Начальные условия
u0 = [7.0, 20.0] # Начальные популяции: x0 = 7, y0 = 20

# Временной интервал
tspan = (0.0, 30.0) # Интервал интегрирования от 0 до 30

# Решение дифференциального уравнения
prob = ODEProblem(pred_prey!, u0, tspan)
sol = solve(prob)

# Построение графика популяций
plot(sol, xlabel="Время", ylabel="Популяция",
     title="Модель Хищник-Жертва")

# Построение фазового портрета
plot(sol[1,:], sol[2,:], xlabel="Численность жертв",
     ylabel="Численность хищников", title="Фазовый Портрет")
savefig("predator_prey_model.png")