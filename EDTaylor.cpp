#include<iostream>
#include<fstream>
#include<cmath>
double y_i;
double analitica(double x){
    double y;
    y=(-1)/(((x*x/4)+(1/2*x)-1));
    return y;
}

double d_f(double x,double y){
return 0.5*(1+x)*pow(y,2);
}

double d_2_f(double x, double y){
return 0.5*pow(y,2)+(1+x)*y*d_f(x,y);
}

double d_3_f(double x,double y){
return 2*y*d_f(x,y)+(1+x)*pow(d_f(x,y),2)+(1+x)*y*d_2_f(x,y);
}

double sol_taylor(double x0,double y0, double h){
    y_i=y0+h*d_f(x0,y0)+0.5*pow(h,2)*d_2_f(x0,y0)+(1/6)*pow(h,3)*d_3_f(x0,y0);
    return y_i;
}

int main(){
    double x_0, y_0;
    std::cout<<"Bienvenido a el metodo EDO de Taylor"<<std::endl;
    std::cout<<"Ingrese el valor de x_0:"<<std::endl;
    std::cin>>x_0;
    std::cout<<"Ingrese el valor de y_0:"<<std::endl;
    std::cin>>y_0;
double xi=x_0;
double y_t=y_0;
    //Recuerde que el intervalo es [0,1]
    //Metodo de Taylor para h=0.01
    //Almacenar en un archivo .txt o .dat junto con la solución exacta.
std::ofstream datafile("resultados.dat");
while(xi<1){
    float exacto=analitica(xi);
    y_t = sol_taylor(xi,y_t,0.01);
    std::cout<<xi<<" "<<y_t<<" "<<exacto<<std::endl;
    datafile<<xi<<" "<<y_t<<" "<<exacto<<std::endl;
    xi=xi+0.01;
    
    
    
    
    
    }
datafile.close();


    //Realizar el metodo de Taylor con h=0.01 y h=0.001




    //El script para graficar. Y exporte la grafica en png o jpg, 
    //hagale un toque personal

    std::ofstream scriptFile1("Grafico_aproximacion_vs_Exacta.gp");
    scriptFile1<<"set term png\n";
    scriptFile1<<"set output 'Grafico_aproximacion_vs_Exacta.png'\n";
    scriptFile1<<"set xlabel 'xi'\n";
    scriptFile1<<"set ylabel 'yi'\n";
    scriptFile1<<"set logscale x\n";
    scriptFile1<<"plot 'resultados.dat'  u 1:3 w l title 'Funcion exacta'\n";
    scriptFile1.close();

system("gnuplot Grafico_aproximacion_vs_Exacta.gp");
double uwu= sol_taylor(0,0.57,0.001);
std::cout<<uwu<<std::endl;
return 0;
}