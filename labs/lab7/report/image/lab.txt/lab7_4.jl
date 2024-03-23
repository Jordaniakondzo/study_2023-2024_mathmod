using Plots
using DifferentialEquations

# Коэффициенты для первого уравнения
a1 = 0.64
a2 = 0.00004

# Коэффициенты для второго уравнения
b1 = 0.00007
b2 = 0.7

# Коэффициенты для третьего уравнения
c1 = 0.4
c2 = 0.3
N = 1403    # Общее количество людей в популяции
n0 = 9      # Начальное количество осведомленных людей

# Определение функции модели ОДУ
function odn_f(du, u, p, t)
    du[1] = (a1 + a2*u[1]) * (N - u[1])    # Модель 1
    du[2] = (b1 + b2*u[1]) * (N - u[1])    # Модель 2
    du[3] = (c1*t + c2*sin(2t)*u[1]) * (N - u[1]) # Модель 3
end

u0 = [n0, n0, n0]  # Начальные условия для каждой модели
tspan = (0.0, 30.0)  # Интервал времени для симуляции

# Создание и решение проблемы ОДУ
prob = ODEProblem(odn_f, u0, tspan)
sol = solve(prob, dtmax=0.1)

# Извлечение данных из решения
N1 = [u[1] for u in sol.u]  # Результаты для первой модели
N2 = [u[2] for u in sol.u]  # Результаты для второй модели
N3 = [u[3] for u in sol.u]  # Результаты для третьей модели
T = [t for t in sol.t]      # Временные метки

# Построение графиков изменения числа осведомленных людей во времени для каждой модели
plt = plot(layout=(1, 3), dpi=300, legend=:bottomright)
plot!(plt[1], T, N1, title="Модель 1", label="Уравнение 1", color=:blue)
plot!(plt[2], T, N2, title="Модель 2", label="Уравнение 2", color=:green)
plot!(plt[3], T, N3, title="Модель 3", label="Уравнение 3", color=:red)

# Сохранение графика
savefig("advertising_models.png")

# Проследите решения для каждой модели
p1 = plot(T, N1, label="Модель 1", color=:blue, xlabel="Время", ylabel="Количество информированных людей", title="Распространение рекламы Модель 1")
p2 = plot(T, N2, label="Модель 2", color=:green, xlabel="Время", ylabel="Количество информированных людей", title="Распространение рекламы Модель 2")
p3 = plot(T, N3, label="Модель 3", color=:red, xlabel="Время", ylabel="Количество информированных людей", title="Распространение рекламы Модель 3")

# Показать графики в новом окне
display(p1)
display(p2)
display(p3)

# Сохранение графиков в файлы
savefig(p1, "model1.png")
savefig(p2, "model2.png")
savefig(p3, "model3.png")

# Функция для расчета производной модели 2
function derivative_model2(n, t, b1, b2, N)
    # Вычисление скорости распространения рекламы в момент времени t
    return (b1 + b2*n)*(N - n)
end

# Вычисление скорости распространения для каждого момента времени в решении модели 2
speeds = [derivative_model2(sol1.u[i][2], sol1.t[i], b1, b2, N) for i in 1:length(sol1.t)]

# Поиск индекса максимальной скорости распространения
max_speed_index = argmax(speeds)

# Время, соответствующее максимальной скорости распространения
max_speed_time = sol1.t[max_speed_index]

println("Время, когда скорость распространения рекламы максимальна, t = $max_speed_time")