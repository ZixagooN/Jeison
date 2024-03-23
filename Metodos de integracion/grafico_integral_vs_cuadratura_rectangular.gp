set term png
set output 'grafico_integral_vs_cuadratura_rectangular.png'
set xlabel 'N'
set ylabel 'Integral'
set logscale x
plot 'resultados_cuadratura_rectangular.dat' u 1:2 w l title 'Integral aproximada', '' u 1:3 w l title 'Integral Simpson', '' u 1:4 w l title 'Integral Trapecio', '' u 1:5 w l title 'Valor exacto'
