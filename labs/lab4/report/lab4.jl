1. Колебания гармонического осциллятора без затуханий и без действий внешней силы

using DifferentialEquations
using Plots

# Определение системы ОДУ для первого случая
# Колебания гармонического осциллятора без затухания и без действия внешней силы
function harmonic_oscillator_1(du, u, p, t)
    x, v = u  # u[1] это положение x, u[2] это скорость v
    du[1] = v             # dx/dt = v
    du[2] = -5*x          # dv/dt = -5x (без затухания, без внешней силы)
end

# Начальные условия и интервал времени
u0 = [0.0, 1.0]          # x0 = 0, v0 = 1
tspan = (0.0, 30.0)      # От t = 0 до t = 30
dt = 0.05                # Шаг времени для сохранения результатов

# Создание проблемы ОДУ
prob1 = ODEProblem(harmonic_oscillator_1, u0, tspan)

# Решение проблемы ОДУ
sol1 = solve(prob1, saveat=dt)

# Построение фазового портрета
plot(sol1, idxs=(1,2), xlabel="Положение (x)", ylabel="Скорость (v)", title="Фазовый портрет гармонического осциллятора 1")

2. Колебания гармонического осциллятора с затуханием и без действий внешней силы

using DifferentialEquations
using Plots

# Определение системы ОДУ для второго случая
# Колебания гармонического осциллятора с затуханием и без действия внешней силы
function harmonic_oscillator_2(du, u, p, t)
    x, v = u  # u[1] это положение x, u[2] это скорость v
    du[1] = v             # dx/dt = v, первое уравнение системы
    du[2] = -2*v - 5*x    # dv/dt = -2v - 5x, второе уравнение системы с затуханием
end

# Начальные условия и интервал времени
u0 = [0.0, 1.0]          # x0 = 0, v0 = 1, начальное положение и начальная скорость
tspan = (0.0, 30.0)      # От t = 0 до t = 30, временной интервал для решения ОДУ
dt = 0.05                # Шаг времени для сохранения результатов

# Создание проблемы ОДУ для решения
prob2 = ODEProblem(harmonic_oscillator_2, u0, tspan)

# Решение проблемы ОДУ
sol2 = solve(prob2, saveat=dt)

# Построение фазового портрета
plot(sol2, idxs=(1,2), xlabel="Положение (x)", ylabel="Скорость (v)", title="Фазовый портрет гармонического осциллятора 2")

3. Колебания гармонического осциллятора с затуханием и с действиями внешней силы :

using DifferentialEquations
using Plots

# Определение системы ОДУ для третьего случая
# Колебания гармонического осциллятора с затуханием и под действием внешней силы
function harmonic_oscillator_3(du, u, p, t)
    x, v = u  # u[1] это положение x, u[2] это скорость v
    du[1] = v                          # dx/dt = v, первое уравнение системы
    du[2] = -4*v - x + sin(14*t)       # dv/dt = -4v - x + sin(14t), второе уравнение системы с затуханием и внешней силой
end

# Начальные условия и интервал времени
u0 = [0.0, 1.0]         # x0 = 0, v0 = 1, начальное положение и начальная скорость
tspan = (0.0, 30.0)     # От t = 0 до t = 30, временной интервал для решения ОДУ
dt = 0.05               # Шаг времени для сохранения результатов

# Создание проблемы ОДУ для решения
prob3 = ODEProblem(harmonic_oscillator_3, u0, tspan)

# Решение проблемы ОДУ
sol3 = solve(prob3, saveat=dt)

# Построение фазового портрета
plot(sol3, idxs=(1,2), xlabel="Положение (x)", ylabel="Скорость (v)", title="Фазовый портрет гармонического осциллятора 3")