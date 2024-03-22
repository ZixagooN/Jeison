#include<iostream>
#include<fstream>
#include<math>

double analitica(double x, double y){
    double y;
    y=(-1)/((x*x*0.25)+(0.5*x)-1)
    y=FORMULA;
    return y;
}

double d_f(double x,double y){
    return(0.5*(1+x)*pow(y,2))
}

double d_2_f(double x,double y){
    return(0.5)*pow(y,2)+(1+x)*y*d_f(x,y)
}

double d_3_f(double x, double y){
    return(2*y*d_f(x,y)+(1+x)*pow(d_f(x,y),2)+(1+x)*y*d_2_f(x,y))
}

double sol_taylor(double x0, double y0, double h){
    double y_i;
    y_i=y0+h*d_f(x0,y0)+0.5*pow(h,2)*d_2_f(x0,y0)+(1/6)*pow(h,3)*d_3_f(x0,y0);  //es y_i+1
    return y_i;
}

int main(){
    std::cout<<"Bienvenido a el metodo EDO de Taylor"<<std::endl;
    std::cout<<"Ingrese el valor de x_0:"<<std::endl;
    std::cin>>x0;
    std::cout<<"Ingrese el valor de y_0:"<<std::endl;
    std::cin>>y0;


//Intervalo de operacion es [0,1]

//realizar metodo de taylor con h=0.01 y h=0.001

//almacenar en un archivo .txt o .dat junto con la solucion exacta

//Script para graficar

return 0
}