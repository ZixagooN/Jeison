//by Jeison Ramos
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
using namespace std;
double v0=0;
double m=70;
double g=9.81;
double rho=1.204;
double A=0.6;
double delta=0.8;
double b=(rho*A*delta)/(2);
double vlim=sqrt((m*g)/b);
double t0=0;
double c=g/vlim;
double tf;
double h=0.01;
double k1, k2, k3, k4;
double opcion;


double dv(double v){
    return b/m*pow(v,2)-g;
}


double d_2_v(double v){
return 2*b*v*dv(v)/m;
}

double d_3_v(double v){
return 2*b*(dv(v)*dv(v)+v*d_2_v(v))/m;
}

double sol_taylor(double v0, double h){
    return v0+h*dv(v0)+0.5*pow(h,2)*d_2_v(v0)+(1/6)*pow(h,3)*d_3_v(v0);;
}


double analitica(double t){
    return vlim*(((v0-vlim)*exp((c)*(t-t0))+(v0+vlim)*exp((-c)*(t-t0)))
    /((v0-vlim)*exp((c)*(t-t0))-(v0+vlim)*exp((-c)*(t-t0))));
}

double ti=t0;
double vi=v0;

//Solucion Analitica
double SolAnalitica(){
std::ofstream datafile("resultados.dat");
for(double j=0; j<=tf; j+=h){
    double integral = analitica(j);
    datafile<<j<<" "<<integral<<std::endl;
    }
datafile.close();
std::ofstream scriptFile1("grafica_aproximada.gp");
    scriptFile1<<"set term png\n";
    scriptFile1<<"set output 'grafica_aproximada.png'\n";
    scriptFile1<<"set xlabel 't'\n";
    scriptFile1<<"set ylabel 'v'\n";
    scriptFile1<<"plot 'resultados.dat' u 1:2 w l lc rgb 'blue' title 'Solución Analitica'\n";
    scriptFile1.close();
return 0;
}



//Grafica de Taylor
double Taylor(){
    std::ofstream datafile("resultados.dat");
        while(ti<tf){
        double vm=-1*vi;
        double solucion = analitica(ti);
        double error = abs(solucion-vm);
        datafile<<ti<<" "<<vm<<" "<<error<<" "<<solucion<<std::endl;
        double v_siguiente = sol_taylor(vi,h);
        vi=v_siguiente;
        ti+=h;
        }
    datafile.close();
        std::ofstream scriptFile1("grafica_aproximada.gp");
        scriptFile1<<"set term png\n";
        scriptFile1<<"set output 'grafica_aproximada.png'\n";
        scriptFile1<<"set xlabel 't'\n";
        scriptFile1<<"set ylabel 'v'\n";
        scriptFile1<<"plot 'resultados.dat' u 1:2 w l lc rgb 'red' title 'Solución Aproximada Taylor', 'resultados.dat' u 1:4 w l lc rgb 'blue' title 'Solución Analitica'\n";
        scriptFile1.close();

        std::ofstream scriptFile2("grafica_aproximada2.gp");
        scriptFile2<<"set term png\n";
        scriptFile2<<"set output 'Error del metodo.png'\n";
        scriptFile2<<"set xlabel 't'\n";
        scriptFile2<<"set ylabel 'v'\n";
        scriptFile2<<"plot 'resultados.dat' u 1:3 w l lc rgb 'blue' title 'Error Taylor'\n";
        scriptFile2.close();
        system("gnuplot grafica_aproximada2.gp");
        return 0;
}

//Grafica de Euler
double Euler(){
    double D=(tf-t0)/h;
std::ofstream datafile("resultados.dat");
    for(int i=0; i<=D; ++i){
        double ti=ti+h;
        double vi= vi +h*dv(vi);
        double vm= -1*vi;
        double solucion = analitica(ti);
        double error = abs(solucion-vm);
        datafile<<ti<<" "<<vm<<" "<<error<<" "<<solucion<<std::endl;
        }
datafile.close();
    std::ofstream scriptFile1("grafica_aproximada.gp");
    scriptFile1<<"set term png\n";
    scriptFile1<<"set output 'grafica_aproximada.png'\n";
    scriptFile1<<"set xlabel 't'\n";
    scriptFile1<<"set ylabel 'v'\n";
    scriptFile1<<"plot 'resultados.dat' u 1:2 w l lc rgb 'red' title 'Solución Aproximada Euler', '' u 1:4 w l lc rgb 'blue' title 'Solución Analitica'\n";
    scriptFile1.close();

    std::ofstream scriptFile2("grafica_aproximada2.gp");
    scriptFile2<<"set term png\n";
    scriptFile2<<"set output 'Error del metodo.png'\n";
    scriptFile2<<"set xlabel 't'\n";
    scriptFile2<<"set ylabel 'v'\n";
    scriptFile2<<"plot 'resultados.dat' u 1:3 w l lc rgb 'blue' title 'Error Euler'\n";
    scriptFile2.close();
    system("gnuplot grafica_aproximada2.gp");
    return 0;
}

//Grafica de Euler Mejorado
double EulerM(){
std::ofstream datafile("resultados.dat");
    while (ti<tf){
    double y_asterisco=vi+h*dv(vi);
    double v_siguiente=vi+0.5*h*(dv(vi)+dv(y_asterisco));
    vi=v_siguiente;
    double vm=-1*vi;
    double solucion = analitica(ti);
    double error = abs(solucion-vm);
    datafile<<ti<<" "<<vm<<" "<<solucion<<" "<<error<<std::endl;
    ti+=h;}
datafile.close();
    std::ofstream scriptFile1("grafica_aproximada.gp");
    scriptFile1<<"set term png\n";
    scriptFile1<<"set output 'grafica_aproximada.png'\n";
    scriptFile1<<"set xlabel 't'\n";
    scriptFile1<<"set ylabel 'v'\n";
    scriptFile1<<"plot 'resultados.dat' u 1:2 w l lc rgb 'red' title 'Solución Aproximada Euler Mejorado', '' u 1:3 w l lc rgb 'blue' title 'Solución Analitica'\n";
    scriptFile1.close();

    std::ofstream scriptFile2("grafica_aproximada2.gp");
    scriptFile2<<"set term png\n";
    scriptFile2<<"set output 'Error del metodo.png'\n";
    scriptFile2<<"set yrange [-0.1 : 0.1]\n";
    scriptFile2<<"set xlabel 't'\n";
    scriptFile2<<"set ylabel 'v'\n";
    scriptFile2<<"plot 'resultados.dat' u 1:4 w l lc rgb 'blue' title 'Error Euler Mejorado'\n";
    scriptFile2.close();
    system("gnuplot grafica_aproximada2.gp");
    return 0;
}

double RungeKutta(){
    std::ofstream datafile("resultados.dat");
        while (ti<tf){
        if (ti+h>tf)
        h = tf-ti;
        k1 = h * dv(vi);
        k2 = h * dv(vi + 0.5 * k1);
        k3 = h * dv(vi + 0.5 * k2);
        k4 = h * dv(vi + k3);
        double vm=-1*(vi);
        double solucion = analitica(ti);
        double error = solucion-vm;
        datafile<<ti<<" "<<vm<<" "<<solucion<<" "<<error<<std::endl;
        vi = vi + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        ti += h;
        }
    datafile.close();
    std::ofstream scriptFile1("grafica_aproximada.gp");
    scriptFile1<<"set term png\n";
    scriptFile1<<"set output 'grafica_aproximada.png'\n";
    scriptFile1<<"set xlabel 't'\n";
    scriptFile1<<"set ylabel 'v'\n";
    scriptFile1<<"plot 'resultados.dat' u 1:2 w l lc rgb 'red' title 'Solución Aproximada RungeKutta', '' u 1:3 w l lc rgb 'blue' title 'Solución Analitica'\n";
    scriptFile1.close();

    std::ofstream scriptFile2("grafica_aproximada2.gp");
    scriptFile2<<"set term png\n";
    scriptFile2<<"set output 'Error del metodo.png'\n";
    scriptFile1<<"set autoscale y\n";
    scriptFile2<<"set xlabel 't'\n";
    scriptFile2<<"set ylabel 'v'\n";
    scriptFile2<<"plot 'resultados.dat' u 1:4 w l lc rgb 'blue' title 'Error RungeKutta'\n";
    scriptFile2.close();
    system("gnuplot grafica_aproximada2.gp");
    return 0;
}


int main(){
std::cout<<"Aproximaciones Numéricas del problema del paracaidista"<<std::endl;
std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
std::cout<<"1. Grafica Analítica."<<std::endl;
std::cout<<"2. Aproximación Numérica Metodo de Taylor, Analítica y Error Absoluto."<<std::endl;
std::cout<<"3. Aproximación Numérica Metodo de Euler, Analítica y Error Absoluto."<<std::endl;
std::cout<<"4. Aproximación Numérica Metodo de Euler Mejorado, Analítica y Error Absoluto."<<std::endl;
std::cout<<"5. Aproximación Numérica Metodo de Runge Kutta, Analítica y Error Absoluto."<<std::endl;
std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
std::cout<<"Ingrese su opción"<<std::endl;
std::cin>>opcion;
std::cout<<"Ingrese el tiempo de simulación"<<std::endl;
std::cin>>tf;
system("clear");
std::cout<<"Ya se ha generado la grafica de la simulación"<<std::endl;

if(opcion==1){
    SolAnalitica();
}
else if (opcion==2)
{
    Taylor();
}
else if (opcion==3)
{
    Euler();
}
else if (opcion==4)
{
    EulerM();
    
}
else if (opcion==5)
{
    RungeKutta();
}
else{
std::cout<<"Opcion Invalida"<<std::endl;
}

system("gnuplot grafica_aproximada.gp");
return 0;
}
