1. Решение для случая 1:

using Plots
using DifferentialEquations

# Параметры модели и начальные условия
p_cr = 47 
N = 50 
q = 1 
tau1 = 33 
tau2 = 27
p1 = 9.7 
p2 = 11.7 
M01 = 7.7
M02 = 9.7

# Вычисление коэффициентов для уравнений
a1 = p_cr/(tau1*tau1*p1*p1*N*q)
a2 = p_cr/(tau2*tau2*p2*p2*N*q)
b = p_cr/(tau1*tau1*tau2*tau2*p1*p1*p2*p2*N*q)
c1 = (p_cr-p1)/(tau1*p1)
c2 = (p_cr-p2)/(tau2*p2)

# Функция дифференциальных уравнений
function ode_f(du, u, p, t)
    # Система ОДУ для двух фирм
    m1, m2 = u
    du[1] = (c1/c1)*u[1]-(a1/c1)*u[1]*u[1]-(b/c1)*u[1]*u[2]
    du[2] = (c2/c1)*u[2]-(a2/c1)*u[2]*u[2]-(b/c1)*u[1]*u[2]
end

# Начальные условия и интервал времени для решения ОДУ
u0 = [M01, M02]
tspan =(0.0, 30.0)
prob1 = ODEProblem(ode_f, u0, tspan)
sol1 = solve(prob1, dtmax=0.1)

# Извлечение данных из решения
M1 = [u[1] for u in sol1.u]
M2 = [u[2] for u in sol1.u]
T = [t for t in sol1.t]

# Построение графиков решений
plt = plot(layout=(1), dpi=300, legend=true)
plot!(plt[1], T, M1, title="График M для случая 1", label="Решение уравнения 1", color=:blue)
plot!(plt[1], T, M2, label="Решение уравнения 2", color=:green)

# Сохранение графика       
savefig("lab8_1.png")

2. Решение для случая 2:

using Plots
using DifferentialEquations

# Задаем параметры модели
p_cr = 47 
N = 50 
q = 1 
tau1 = 33 
tau2 = 27
p1 = 9.7 
p2 = 11.7 
M01 = 7.7
M02 = 9.7

# Вычисляем коэффициенты для дифференциальных уравнений
a1 = p_cr/(tau1*tau1*p1*p1*N*q)
a2 = p_cr/(tau2*tau2*p2*p2*N*q)
b = p_cr/(tau1*tau1*tau2*tau2*p1*p1*p2*p2*N*q)
c1 = (p_cr-p1)/(tau1*p1)
c2 = (p_cr-p2)/(tau2*p2)

# Функция для системы ОДУ
function ode_f(du, u, p, t)
    # Система ОДУ для двух фирм
    m1, m2 = u
    du[1] = (c1/c1)*u[1]-(a1/c1)*u[1]*u[1]-((b/c1)+0.00044)*u[1]*u[2]
    du[2] = (c2/c1)*u[2]-(a2/c1)*u[2]*u[2]-(b/c1)*u[1]*u[2]
end

# Начальные условия и решение системы ОДУ
u0 = [M01, M02]
tspan =(0.0, 30.0)
prob1 = ODEProblem(ode_f, u0, tspan)
sol1 = solve(prob1, dtmax=0.1)

# Извлекаем решения для M1 и M2
M1 = [u[1] for u in sol1.u]
M2 = [u[2] for u in sol1.u]
T = [t for t in sol1.t]

# Построение графиков изменений M1 и M2
plt = plot(T, M1, label="Решение уравнения 1", color=:red, title="График M1", xlabel="Время", ylabel="M1")
plot!(T, M2, label="Решение уравнения 2", color=:blue, title="График M2", xlabel="Время", ylabel="M2")

# Сохраняем графики в файл
savefig("lab8_2.png")