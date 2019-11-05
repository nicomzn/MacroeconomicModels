###############################################################################
#       Simulacion del modelo de ciclos de Goodwin usando Julia Language
#                       Catedra de Macroeconomia II - UNLP
#                                 Nicolas Monzon
#
# INTRODUCCION
#
# El siguiente codigo es utilizado para realizar las simulacion numerica y
# grafica del modelo de ciclos de Goodwin (1967). Para realizar las simulaciones
# se utilizaran las estimaciones realizadas por Harvie (2000) que utiliza datos
# para los paises miembros de la OECD para el periodo 1951-1994
#
#
# REFERENCIAS
#
# Harvie, D. (2000). Testing Goodwin: growth cycles in ten OECD countries.
# Cambridge Journal of Economics, 24(3), 349-376.
#
###############################################################################
# Directorio de trabajo
cd("D:/UNLP - Macro II/2019/Goodwin/Julia/Resultados")

# Cargo el paquete utilizado para resolver ODEs
using DifferentialEquations

# Armo la funcion que resuelve las formas semi-reducidas del modelo de goodwin
function goodwin(dy,y,parameter,t)
  v, u = y
  α, β, σ, γ, ρ = parameter
  dy[1] = dv = (1/σ - (α + β))*v - (1/σ)*u*v
  dy[2] = du = -(α + γ)*u + ρ*v*u
end

# Para resolver el modelo es necesario darle valores a los parametros declarados
# y establecer el periodo de tiempo que vamos a utilizar

# Para calibrar parametros y el tiempo del ciclo utilizamos el paper de Harvie (2000)
# Vamos a tener un vector de parametros y una duracion del ciclo diferente para
# cada pais => T = 2*π/((α+γ)*(1/σ-(α+β)))^(1/2)

# Australia
#          α       Β       σ       γ      ρ
p_aus = [0.0166, 0.0226, 2.4994, 62.36, 67.1]
T_aus = 2*π /((0.0166+62.36)*(1/2.4994-(0.0166+0.0226)))^(1/2)
tspan_aus = (0.0,T_aus)
#          v     u
u0_aus = [0.93, 0.9]

# Canada
#          α       Β       σ       γ      ρ
p_can = [0.0160, 0.0259, 1.5698, 59.01, 65.32]
T_can = 2*π /((0.0160+59.01)*(1/1.5698-(0.0160+0.0259)))^(1/2)
tspan_can = (0.0,T_can)
#          v     u
u0_can = [0.90, 0.93]

# Finlandia
#          α       Β       σ       γ      ρ
p_fin = [0.0303, 0.008, 3.1396, 32.00, 36.57]
T_fin = 2*π /((0.0303+32.00)*(1/3.1396-(0.0303+0.008)))^(1/2)
tspan_fin = (0.0,T_fin)
#          v     u
u0_fin = [0.88, 0.89]

# Francia
#          α       Β       σ       γ      ρ
p_fra = [0.0364, 0.0076, 1.7974, 54.85, 62.01]
T_fra = 2*π /((0.0364+54.85)*(1/1.79746-(0.0364+0.0076)))^(1/2)
tspan_fra = (0.0,T_fra)
#          v     u
u0_fra = [0.89, 0.92]

# Alemania
#          α       Β       σ       γ      ρ
p_ale = [0.0329, 0.0041, 2.4941, 85.49, 65.55]
T_ale = 2*π /((0.0329+85.49)*(1/2.4941-(0.0329+0.0041)))^(1/2)
tspan_ale = (0.0,T_ale)
#          v     u
u0_ale = [1.13, 0.90]

# Grecia
#          α       Β       σ       γ      ρ
p_gre = [0.0401, 0.0035, 3.0292, 46.02, 53.48]
T_gre = 2*π /((0.0401+46.02)*(1/3.0292-(0.0401+0.0035)))^(1/2)
tspan_gre = (0.0,T_gre)
#          v     u
u0_gre = [0.86, 0.87]

# Italia
#          α       Β       σ       γ      ρ
p_ita = [0.0460, 0.0049, 3.3527, 71.24, 81.97]
T_ita = 2*π /((0.0460+71.24)*(1/3.3527-(0.0460+0.0049)))^(1/2)
tspan_ita = (0.0,T_ita)
#          v     u
u0_ita = [0.87, 0.83]

# Noruega
#          α       Β       σ       γ      ρ
p_nor =[0.0262, 0.0134, 3.6710, 118.07, 122.43]
T_nor = 2*π /((0.0262+118.07)*(1/3.6710-(0.0262+0.0134)))^(1/2)
tspan_nor = (0.0,T_nor)
#          v     u
u0_nor = [0.96, 0.86]

# UK
#          α       Β       σ       γ      ρ
p_uk = [0.0221, 0.0036, 2.5694, 18.54, 21.90]
T_uk = 2*π /((0.0221+18.54)*(1/2.5694-( 0.0221+0.0036)))^(1/2)
tspan_uk = (0.0,T_uk)
#          v     u
u0_uk = [0.85, 0.93]

# Defino el problema para cada pais
prob_aus = ODEProblem(goodwin,u0_aus,tspan_aus,p_can)
prob_can = ODEProblem(goodwin,u0_can,tspan_can,p_can)
prob_fin = ODEProblem(goodwin,u0_fin,tspan_fin,p_fin)
prob_fra = ODEProblem(goodwin,u0_fra,tspan_fra,p_fra)
prob_ale = ODEProblem(goodwin,u0_ale,tspan_ale,p_ale)
prob_gre = ODEProblem(goodwin,u0_gre,tspan_gre,p_gre)
prob_ita = ODEProblem(goodwin,u0_ita,tspan_ita,p_ita)
prob_nor = ODEProblem(goodwin,u0_nor,tspan_nor,p_nor)
prob_uk  = ODEProblem(goodwin,u0_uk,tspan_uk,p_uk)

# Defino la solucion para cada problema
gw_aus = solve(prob_aus,Tsit5(),saveat=0.1)
gw_can  = solve(prob_can ,Tsit5(),saveat=0.1)
gw_fin = solve(prob_fin,Tsit5(),saveat=0.1)
gw_fra = solve(prob_fra,Tsit5(),saveat=0.1)
gw_ale = solve(prob_ale,Tsit5(),saveat=0.1)
gw_gre = solve(prob_gre,Tsit5(),saveat=0.1)
gw_ita = solve(prob_ita,Tsit5(),saveat=0.1)
gw_nor = solve(prob_nor,Tsit5(),saveat=0.1)
gw_uk = solve(prob_uk,Tsit5(),saveat=0.1)

# Cargo el paquete de graficos
using Plots

# Tendencias de cada solucion
trends_aus = plot(gw_aus, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_can = plot(gw_can, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_fin = plot(gw_fin, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_fra = plot(gw_fra, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_ale = plot(gw_ale, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_gre = plot(gw_gre, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_ita = plot(gw_ita, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_nor = plot(gw_nor, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")
trends_uk  = plot(gw_uk, title="Tendencias del ciclo",
            label=["Tasa de empleo" "Part de los salarios en el producto"],
            xlabel = "Duracion del Ciclo")

# Diagramas de fase
phase_aus = plot(gw_aus, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_can = plot(gw_can, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_fin = plot(gw_fin, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_fra = plot(gw_fra, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_ale = plot(gw_ale, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_gre = plot(gw_gre, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_ita = plot(gw_ita, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_nor = plot(gw_nor, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])
phase_uk  = plot(gw_uk, plotdensity=10000, vars=(1,2),
title="Diagrama de Fases", label=["Ciclo eliptico de u y v"])

# Combino ambos graficos y los exporto
aus = plot(trends_aus, phase_aus)
png(aus, "Australia")
can = plot(trends_can, phase_can)
png(can, "Canada")
fin = plot(trends_fin, phase_fin)
png(fin, "Finlandia")
fra = plot(trends_fra, phase_fra)
png(fra, "Francia")
ale = plot(trends_ale, phase_ale)
png(ale, "Alemania")
gre = plot(trends_gre, phase_gre)
png(gre, "Grecia")
ita = plot(trends_ita, phase_ita)
png(ita, "Italia")
nor = plot(trends_nor, phase_nor)
png(nor, "Noruega")
uk  = plot(trends_uk, phase_uk)
png(uk, "UK")


# Animacion del Diagrama de fase
gif_phase_aus = animate(gw_aus,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_can = animate(gw_can,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_fin = animate(gw_fin,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_fra = animate(gw_fra,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_ale = animate(gw_ale,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_gre = animate(gw_gre,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_ita = animate(gw_ita,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_nor = animate(gw_nor,lw=3 ,every=1, plotdensity=10000, vars=(1,2))
gif_phase_uk  = animate(gw_uk,lw=3 ,every=1, plotdensity=10000, vars=(1,2))


###############################################################################
