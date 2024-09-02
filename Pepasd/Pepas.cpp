#include<iostream>
#include<fstream>


int t=0;
int x1;
int x2;
int x3;
int x4;
int posicion1(double t){
    return 4*t+7;
}

int posicion2(double t){
    return 9*t+5;
}


int main(){
std::ofstream datafile("resultados.dat");
while(t<10){
    x1=posicion1(t);
    x2=posicion2(t);
    x3=posicion1(5*t);
    x4=posicion2(5*t);
    datafile<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<std::endl;
    t=t+1;

}

datafile.close();

std::ofstream scriptFile1("x1_vs_x2.gp");
scriptFile1<<"set term pdf\n";
scriptFile1<<"set title 'x_{1} vs x_{2} de t [0-9]'\n";
scriptFile1<<"set output 'grafica x1 vs x2.pdf'\n";
scriptFile1<<"set xlabel 'x_{1}'\n";
scriptFile1<<"set ylabel 'x_{2}'\n";
scriptFile1<<"plot 'resultados.dat' u 1:2 w p pt 7 linecolor rgb 'blue' title 'x1 vs x2' at 5.0, 5.0\n";
scriptFile1.close();

std::ofstream scriptFile2("x1_vs_x2_linea.gp");
scriptFile2<<"set term pdf\n";
scriptFile2<<"set title 'x_{1} vs x_{2} de t [0-9]'\n";
scriptFile2<<"set output 'grafica linea x1 vs x2.pdf'\n";
scriptFile2<<"set xlabel 'x_{1}'\n";
scriptFile2<<"set ylabel 'x_{2}'\n";
scriptFile2<<"plot 'resultados.dat' u 1:2 w lp pt 7 linecolor rgb 'purple' lw 3 title 'x1 vs x2' at 5.0, 5.0\n";
scriptFile2.close();

std::ofstream scriptFile3("x1_vs_x2_50.gp");
scriptFile3<<"set term pdf\n";
scriptFile3<<"set title 'x_{1} vs x_{2} de t [0-50]'\n";
scriptFile3<<"set output 'grafica 50 seg x1 vs x2.pdf'\n";
scriptFile3<<"set xlabel 'x_{1}'\n";
scriptFile3<<"set ylabel 'x_{2}'\n";
scriptFile3<<"plot 'resultados.dat' u 3:4 w p pt 7 linecolor rgb 'blue' title 'x1 vs x2' at 5.0, 5.0\n";
scriptFile3.close();

std::ofstream scriptFile4("x1_vs_x2_linea50.gp");
scriptFile4<<"set term pdf\n";
scriptFile4<<"set title 'x_{1} vs x_{2} de t [0-50]'\n";
scriptFile4<<"set output 'grafica 50 seg linea x1 vs x2.pdf'\n";
scriptFile4<<"set xlabel 'x_{1}'\n";
scriptFile4<<"set ylabel 'x_{2}'\n";
scriptFile4<<"plot 'resultados.dat' u 3:4 w lp pt 7 linecolor rgb 'purple' lw 3 title 'x1 vs x2' at 5.0, 5.0\n";
scriptFile4.close();

system("gnuplot x1_vs_x2.gp");
system("gnuplot x1_vs_x2_linea.gp");
system("gnuplot x1_vs_x2_50.gp");
system("gnuplot x1_vs_x2_linea50.gp");
return 0;
}
