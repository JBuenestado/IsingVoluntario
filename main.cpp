#include <iostream>
#include <thread>
#include "pthread.h"
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <mutex>
#include <fstream>
//  g++ main.cpp -pthread -std=c++17 -o I

using namespace std;

const int MAX = 100000; // para el random
const bool mostrar = false;
const bool randomStart = false;
const bool threadpertemperature = false; // crear un thread en paralelo para ejecutar para cada temperatura. si quieres ejecutarlas
// secuencialmente; false y numeroTemperaturas != 1
const int PMC = 100000; // pasos montecarlo
const int N = 128;       // tamaño matrix
const int numeroTemperaturas = 1;
double temp[numeroTemperaturas];
double tk = 1.5;  // temperatura del sistema, para utilizar este dato, numeroTemperaturas = 1;
int total = 1; // cuantas veces quieres ejecutarlo, multithread, usar con una unica temperatura
ofstream correlacion, datosout, print, printT[numeroTemperaturas];

void calcular(int runs, int &changes, int thread_number, double tk_value) {    
    unsigned int seed = time(NULL) * thread_number;
    cout << "thread: " << thread_number << " . tk: " << tk_value;
    cout << "tk: " << tk << endl;
    tk = tk_value;

    double promedioMagnetizacion, promedioEnergia, calor;
    int i, j, k, n, l;
    int x0, y0, x1, x2, y1, y2;
    int check1, check2, test;
    double autoc[N][N][N];
    double red[N][N], energia, p, xhi, copia[N][N], autocorrelacion[N];
    double calculoMagnetizacion, calculoEnergia;
    test = 0;

    for (int k = 0; k < runs; k++) {
        if (randomStart) {
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    red[i][j] = pow(-1, rand_r(&seed) % (MAX));
                    copia[i][j] = red[i][j];
                }
            }
        } else {
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    red[i][j] = 1;
                    copia[i][j] = red[i][j];
                }
            }
        }

        promedioMagnetizacion = 0;
        calculoEnergia = 0;

        for (i = 0; i < N; i++) {
            autocorrelacion[i] = 0;
            for (j = 0; j < N; j++) {
                for (l = 0; l < N; l++) {
                    autoc[i][j][l] = 0;
                }
            }
        }

        for (int cont = 0; cont < PMC; cont++) {
            for (int w = 0; w < N*N ; w++) {
                x0 = rand_r(&seed) % N;
                y0 = rand_r(&seed) % N;

                x1 = x0;
                y1 = y0;
                x2 = x0;
                y2 = y0;

                if (x0 == N - 1) { //condiciones de contorno
                    x1 = -1;
                }
                if (x0 == 0) {
                    x2 = N;
                }
                if (y0 == N - 1) {
                    y1 = -1;
                }
                if (y0 == 0) {
                    y2 = N;
                }

                energia = 2.0 * red[x0][y0] * (red[x1 + 1][y0] + red[x2 - 1][y0] + red[x0][y1 + 1] + red[x0][y2 - 1]);
                p = min(1.0, exp(-energia / tk)); // calculo de probabilidad
                xhi = (rand_r(&seed) % (MAX));   // MAX;
                if (p >= xhi / MAX) {            // cambio de signo
                    red[x0][y0] = -red[x0][y0];
                }

                if ((w == 0) && (cont % 100 == 0)) {    // calculo cosas voluntario
                    test += 1;
                    calculoMagnetizacion = 0;
                    calculoEnergia += energia;
                    for (i = 0; i < N; i++) {

                        for (j = 0; j < N; j++) {
                            for (l = 0; l < N; l++) {
                                x1 = i + j;
                                if (i + j >= N) { //contorno
                                    x1 = x1 - N;
                                }
                                autoc[i][j][l] += red[j][l] * red[x1][l];
                            }
                            calculoMagnetizacion += red[i][j]; //magnetizacion
                        }
                    }
                    promedioMagnetizacion += abs(calculoMagnetizacion) / (N * N);
                }  //end bucle calculo cosas voluntario
            }
        }     //end bucle iteraciones

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (l = 0; l < N; l++) {
                    autoc[i][j][l] = autoc[i][j][l] / (PMC / 100);
                }
            }
        }

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (l = 0; l < N; l++) {
                    autocorrelacion[i] += autoc[i][j][l] / (N * N);
                }
            }
        }

        promedioMagnetizacion = promedioMagnetizacion / (PMC / 100);
        promedioEnergia = calculoEnergia / (2 * N * (PMC / 100));
        calor = (((pow(calculoEnergia, 2)) / (PMC / 100)) - pow(calculoEnergia / (PMC / 100), 2)) / (N * N * tk);

        if (total == 1) {
            correlacion << "#Valores de autocorrelacion para matriz de tamaño: " << N << " y temperatura: " << tk
                        << " ºK" << endl;
            for (l = 0; l < N; l++) {
                correlacion << l << " " << autocorrelacion[l] << endl;
                print << l << " " << autocorrelacion[l] << endl;
                if(threadpertemperature){
                    printT[thread_number] << l << " " << autocorrelacion[l] << endl;
                }
            }

            correlacion << endl;
            datosout << "Matriz de tamaño: " << N << " y temperatura: " << tk << " ºK" << endl;
            datosout << "Magnetización promedio: " << promedioMagnetizacion << endl;
            datosout << "Energia media: " << promedioEnergia << endl;
            datosout << "Calor Específico: " << calor << " J" << endl << endl;
        }

        if ((total == 1) && (mostrar)) //mostrar las matrices iniciales/finales
        {
            cout << endl
                 << "Red inicial                                            Red final" << endl;
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if (copia[i][j] == 1) {
                        cout << "+";
                    }
                    cout << copia[i][j] << " ";
                }
                cout << "                 ";
                for (j = 0; j < N; j++) {
                    if (red[i][j] == 1) {
                        cout << "+";
                    }
                    cout << red[i][j] << " ";
                }
                cout << endl;
            }

            cout << endl
                 << "Valores Cambiados " << endl;
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if (red[i][j] + copia[i][j] == 2) {
                        cout << "+";
                    }
                    if (red[i][j] + copia[i][j] == 0) {
                        cout << " ";
                    }
                    cout << (red[i][j] + copia[i][j]) / 2 << " ";
                }
                cout << endl;
            }
            cout << endl;
        }

        check1 = 0;
        check2 = 0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (red[i][j] != 1) {
                    check1 = -1;
                }
                if (red[i][j] != -1) {
                    check2 = 1;
                }
            }
        }
        if (check1 == 0 || check2 == 0) {
            changes += 1;
        }

    } //end bucle runs
}

int create_threads(int total_threads) {
    int thread_step = total / total_threads;
    int return_values[total_threads];

    cout << "-------------------Start of Program------------------" << endl;
    if (randomStart) {
        cout << "Matriz inicial aleatoria" << endl;
    } else {
        cout << "Matriz inicial uniforme" << endl;
    }

    cout << endl << "Veces ejecutado: " << total << endl;
    cout << "Numero de Temperaturas para los cuales se ejecuta: " << numeroTemperaturas << endl;
    cout << "Dimension de la Matriz Cuadrada: " << N << endl;
    cout << "Pasos de Montecarlo: " << PMC << endl << endl;

    correlacion.open("autocorrelation.txt");
    datosout.open("datos.txt");

    if(threadpertemperature){ 
        total_threads = numeroTemperaturas;
        for (int i = 0; i < total_threads; ++i) {
            return_values[i] = 0;
        }
        std::chrono::steady_clock::time_point begin_chrono = std::chrono::steady_clock::now();
        cout << "-------------------Start of Iteration------------------" << endl;
        std::vector<std::thread> threads;
        
        for (int i = 0; i < total_threads; ++i){
            if (numeroTemperaturas != 1) {
                tk = temp[i];
            }
            printT[i].open("autocorrelacion" + to_string(tk) + ".txt");

            threads.emplace_back(thread(calcular, thread_step, ref(return_values[i]), i, tk));
        }
            for (auto &thread: threads) {
                thread.join();
            }

            int result = 0;
            for (int i = 0; i < total_threads; ++i) {
                result += return_values[i];
            }
        
        for (int i = 0; i < total_threads; ++i){
            printT[i].close();
        } 
        
        chrono::steady_clock::time_point end_chrono = chrono::steady_clock::now();
            cout << "Time difference All executions (sec) = "
                << (chrono::duration_cast<chrono::microseconds>(end_chrono - begin_chrono).count()) / 1000000.0 << endl;
            cout << "-----------------------------------------------------" << endl;

    } else {

        std::chrono::steady_clock::time_point chronototal = std::chrono::steady_clock::now();
        for (int l = 0; l < numeroTemperaturas; ++l) {
            print.open("autocorrelacion" + to_string(tk) + ".txt");
            if (numeroTemperaturas != 1) {
                tk = temp[l];
            }

            for (int i = 0; i < total_threads; ++i) {
                return_values[i] = 0;
            }

            // cout current time and date
            std::chrono::steady_clock::time_point begin_chrono = std::chrono::steady_clock::now();

            cout << "-------------------Start of Iteration------------------" << endl;
            std::vector<std::thread> threads;

            for (int i = 0; i < total_threads; ++i) { //nuevo
                threads.emplace_back(thread(calcular, thread_step, ref(return_values[i]), i, tk));
            }

            for (auto &thread: threads) {
                thread.join();
            }

            int result = 0;

            for (int i = 0; i < total_threads; ++i) {
                result += return_values[i];
            }

            cout << "Temperatura de los sistemas: " << tk << endl;
            cout << "Combinaciones mezcladas:  " << total - result << endl;
            cout << "Combinaciones uniformes: " << result << endl << endl;

            cout << "-----------------------------------------------------" << endl;
            chrono::steady_clock::time_point end_chrono = chrono::steady_clock::now();
            cout << "Time difference (sec) = "
                << (chrono::duration_cast<chrono::microseconds>(end_chrono - begin_chrono).count()) / 1000000.0 << endl;
            cout << "-----------------------------------------------------" << endl;
            print.close();
        }
        chrono::steady_clock::time_point total_chrono = chrono::steady_clock::now();
         cout << "Time difference All (sec) = "
                << (chrono::duration_cast<chrono::microseconds>(total_chrono - chronototal).count()) / 1000000.0 << endl;
            cout << "-----------------------------------------------------" << endl;

    }
    datosout.close();
    correlacion.close();

    cout << "--------------------End of Program-------------------" << endl;

    return 0;
} // end create threads

int main(int argc, char *argv[]) {
    cout << "=====================================================" << endl;
    int total_threads = 8;
    ifstream datos;
    if (argc > 1) {
        total_threads = atoi(argv[1]);
    }
    if (total_threads > total) {
        total_threads = total;
    }
    if (total % total_threads != 0) { //para que cada hebra tenga el mismo número de ejecuciones al dividir el step
        total += total_threads - total % total_threads;
    }
    datos.open("tk.txt");

    for (int i = 0; i < numeroTemperaturas; i++) {
        datos >> temp[i];
    }
    datos.close();

    cout << "threads: " << total_threads << "\n";
    cout << "runs per thread: " << total / total_threads << endl;
    create_threads(total_threads);

    return 0;
}