#include<iostream>
#include<cmath>
#include<fstream>


//COMPARACION DE METODOS DE INTEGRACION
//BY JEISON RAMOS.

//Función para integrar
double funcion(double x){
    return 100*pow(x,2)*cos(20*x);
}

//Calcuilo de la integral aproximada bajo cuadratura rectangular
double cuadraturaRectandular(double a, double b, int N){
    double h=(b-a)/N;
    double suma=0.0;
    for(int i=0;i<N;i++){
        double xi=a+i*h;
        suma+=funcion(xi);
    }
    return h*suma;
}

//Calcuilo de la integral aproximada metodo de simpson
double Simpson(double a, double b, int N){
   double h=(b-a)/N;
    double suma=0.0;
    double suma1=0.0;
    double x0=0.0;
    for(int i=0;i<N;i++){
        double xi=a+i*h;
                if(i>0 and i%2==0){
            suma+=2*funcion(xi);
        }
        if(i>0 and i%2!=0){
            suma1+=4*funcion(xi);
        }
    }
    return (h/3)*(suma+suma1+funcion(a)+funcion(b));
}

//Calcuilo de la integral aproximada metodo de trapecio
double Trapecio(double a, double b, int N){
   double h=(b-a)/N;
    double suma=0.0;
    double suma1=0.0;
    double x0=0.0;
    for(int i=0;i<N;i++){
        double xi=a+i*h;
                if(i>0){
            suma+=2*funcion(xi);
    }
    }
    return (h/2.0)*(suma+funcion(a)+funcion(b));
}



int main(){

std::ofstream datafile("resultados_cuadratura_rectangular.dat");

for(int N=1; N<=10000; ++N){
    double integral = cuadraturaRectandular(0.0,1.0,N);
    double integral_simpson=Simpson(0.0,1.0,N);
    double integral_exacta=4.7459;
    double integral_trapecio=Trapecio(0.0,1.0,N);
    datafile<<N<<" "<<integral<<" "<<integral_simpson<<" "<<integral_trapecio<<" "<<integral_exacta<<std::endl;
    }
datafile.close();


datafile.close();

std::ofstream scriptFile1("grafico_integral_vs_cuadratura_rectangular.gp");
scriptFile1<<"set term png\n";
scriptFile1<<"set output 'grafico_integral_vs_cuadratura_rectangular.png'\n";
scriptFile1<<"set xlabel 'N'\n";
scriptFile1<<"set ylabel 'Integral'\n";
scriptFile1<<"set logscale x\n";
scriptFile1<<"plot 'resultados_cuadratura_rectangular.dat' u 1:2 w l title 'Integral aproximada', '' u 1:3 w l title 'Integral Simpson', '' u 1:4 w l title 'Integral Trapecio', '' u 1:5 w l title 'Valor exacto'\n";
scriptFile1.close();

//Generar datos para graficas y pintar los rectangulos

std::ofstream dataFile2("datos_cuadratura_rectangular.dat");
int N_plot=10;
double h_plot=(1.0-0.0)/N_plot;
for(double x=0.0; x<=1.0; x+=h_plot){
    dataFile2<<x<<" "<<funcion(x)<<std::endl;
    dataFile2<<x<<" "<<0.0<<std::endl;
    dataFile2<<std::endl;
}
dataFile2.close();

std::ofstream scriptFile2("grafico_funcion_rectangulos_cuadratura_rectangular.gp");
scriptFile2<<"set term png\n";
scriptFile2<<"set output 'grafico_funcion_rectangulos_cuadratura_rectangular.png'\n";
scriptFile2<<"set xlabel 'x'\n";
scriptFile2<<"set ylabel 'f(x)'\n";
scriptFile2<<"set logscale x\n";
scriptFile2<<"plot 'datos_cuadratura_rectangular.dat' w l title 'Funcion', '' w boxes title 'Rectangulos'\n";
scriptFile2.close();
//Ejecutar Gnuplot
system("gnuplot grafico_integral_vs_cuadratura_rectangular.gp");
system("gnuplot grafico_funcion_rectangulos_cuadratura_rectangular.gp");

return 0;
}