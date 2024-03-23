set term png
set output 'grafico_funcion_rectangulos_cuadratura_rectangular.png'
set xlabel 'x'
set ylabel 'f(x)'
set logscale x
plot 'datos_cuadratura_rectangular.dat' w l title 'Funcion', '' w boxes title 'Rectangulos'
