#include <stdio.h>
#include <stdlib.h>
#include <math.h>           // Para funciones matemáticas como tanh()
#include <gsl/gsl_math.h>    
#include <gsl/gsl_vector.h> 

/**********************************************************************/
/*                  MANEJO DE ARCHIVOS                              */
/**********************************************************************/
/*
   Funciones para el manejo de archivos:
   - writeToFile: Abre un archivo en modo escritura.
   - openfile: Define el nombre del archivo (basado en el tiempo) y lo abre.
   - closefile: Cierra el archivo abierto.
*/
void writeToFile(double);
void openfile(double);
void closefile();
char filename[50];  // Nombre del archivo
FILE *file;         // Puntero para manejar el archivo

/**********************************************************************/
/*              PARÁMETROS DE LA SIMULACIÓN                         */
/**********************************************************************/
/*
   Variables globales de simulación:
   - x_final, x_init: Límites espacial del dominio.
   - NX: Número de puntos espaciales.
   - DX: Tamaño del paso en x.
   - t_final, t_init: Límites temporales.
   - NT: Número de pasos de tiempo.
   - DT: Tamaño del paso temporal.
   - kappa: Coeficiente adiabático.
*/
double x_final, x_init, NX, DX;
double t_final, t_init, NT, DT, kappa;

/*
   Variables auxiliares:
   - t: Tiempo actual de la simulación.
   - x: Posición espacial actual.
   - u1tempo: Valor temporal para la presión.
   - u2tempo: Valor temporal para la densidad.
   - u3tempo: Valor temporal para la velocidad.
   - u1i, u2i, u3i: Condiciones en la frontera (x = 0) para mantener el valor fijo.
   - funf1, funf2, funf3: Variables para almacenar los flujos en cada celda.
   - q1, q2, q3: Variables para almacenar las cantidades conservadas (densidad, momento y energía).
*/
double t, x, u2tempo, u3tempo, u1tempo;
double u1i, u2i, u3i;
double funf1, funf2, funf3;
double q1, q2, q3;

/**********************************************************************/
/*           INICIALIZACIÓN DE LOS PARÁMETROS DE SIMULACIÓN           */
/**********************************************************************/
/*
   Función initialize_variables():
   Inicializa los parámetros globales de la simulación.
*/
void initialize_variables() {
    x_final = 10.0;       // Límite derecho del dominio espacial
    x_init = 1.0;         // Límite izquierdo del dominio espacial
    NX = 1000;            // Número total de puntos en el eje x
    DX = (x_final - x_init) / NX; // Tamaño del paso espacial

    t_final = 1.0;        // Tiempo final de la simulación
    t_init = 0.0;         // Tiempo inicial
    NT = 1000;            // Número total de pasos de tiempo
    DT = (t_final - t_init) / NT; // Tamaño del paso temporal

    kappa = 1.4;          // Coeficiente adiabático (valor típico para algunos gases)
}

/**********************************************************************/
/*           INICIALIZACIÓN DE CONDICIONES INICIALES                  */
/**********************************************************************/
/*
   Función initialize():
   Inicializa la malla espacial y establece las condiciones iniciales.
   Se implementa una transición suave entre dos estados usando tanh().
   
   Parámetros:
   - u1in, u2in, u3in: Vectores de GSL que almacenarán la presión, densidad y velocidad
                        en cada punto espacial.
   - tin: Tiempo actual. Si tin == 0 se establecen las condiciones iniciales.
   
   Condiciones:
   - Para x < x0: estado izquierdo (u1_left = 5.0, u2_left = 4.0, u3_left = 5.75).
   - Para x >= x0: estado derecho (u1_right = 0.5, u2_right = 0.5, u3_right = 0.75).
   - La transición se suaviza usando:
         u = u_left + 0.5*(u_right - u_left)*(1 + tanh((x - x0)/epsilon))
   donde:
         x0: centro del dominio.
         epsilon: parámetro que controla el ancho de la transición.
*/
void initialize(gsl_vector *u1in, gsl_vector *u2in, gsl_vector *u3in, double tin) {
    if (tin == 0) {
        openfile(tin);
        
        // Definir centro de la transición y ancho de la misma
        double x0 = (x_init + x_final) / 2.0;
        double epsilon = (x_final - x_init) / 50.0; // Ajustable según necesidad

        for (int i = 0; i < NX; i++) {
            x = x_init + DX * i;
            
            // Se suaviza la transición con tanh()
            u1tempo = 5.0 + 0.5 * (0.5 - 5.0) * (1 + tanh((x - x0) / epsilon));
            u2tempo = 4.0 + 0.5 * (0.5 - 4.0) * (1 + tanh((x - x0) / epsilon));
            u3tempo = 5.75 + 0.5 * (0.75 - 5.75) * (1 + tanh((x - x0) / epsilon));

            gsl_vector_set(u1in, i, u1tempo);
            gsl_vector_set(u2in, i, u2tempo);
            gsl_vector_set(u3in, i, u3tempo);

            // Se escribe la condición inicial en el archivo
            fprintf(file, "%g\t%g\t%g\t%g\n", x,
                    gsl_vector_get(u1in, i),
                    gsl_vector_get(u2in, i),
                    gsl_vector_get(u3in, i));

            // Se guarda el valor en x = 0 para la condición de frontera
            if (i == 0) {
                u1i = gsl_vector_get(u1in, 0);
                u2i = gsl_vector_get(u2in, 0);
                u3i = gsl_vector_get(u3in, 0);
            }
        }
        closefile();
    } else {
        // Para t > 0, se mantiene la condición fija en la frontera x=0
        gsl_vector_set(u1in, 0, u1i);
        gsl_vector_set(u2in, 0, u2i);
        gsl_vector_set(u3in, 0, u3i);
    }
}

/**********************************************************************/
/*                CÁLCULO DE LOS FLUJOS                             */
/**********************************************************************/
/*
   Función fluxes():
   Calcula los flujos en cada celda espacial a partir de las variables primitivas.
   Se asume:
     - u1: presión.
     - u2: densidad.
     - u3: velocidad.
     
   Los flujos se definen como:
     f1 = u2 * u3
     f2 = u2 * u3^2 + u1
     f3 = u2 * u3 * (0.5 * u3^2 + (kappa/(kappa - 1)) * u1/u2)
*/
void fluxes(double funu1, double funu2, double funu3) {
    funf1 = funu2 * funu3;
    funf2 = funu2 * funu3 * funu3 + funu1;
    funf3 = funu2 * funu3 * (0.5 * funu3 * funu3 + (kappa / (kappa - 1)) * funu1 / funu2);
}

/**********************************************************************/
/*         CÁLCULO DE LAS CANTIDADES CONSERVADAS                      */
/**********************************************************************/
/*
   Función charges():
   Calcula las cantidades conservadas a partir de las variables primitivas.
   
   Se definen:
     q1 = densidad (u2)
     q2 = momento = densidad * velocidad (u2 * u3)
     q3 = energía total = 0.5 * u2 * u3^2 + u1/(kappa - 1)
*/
void charges(double funq1, double funq2, double funq3) {
    q1 = funq2;
    q2 = funq2 * funq3;
    q3 = 0.5 * funq2 * funq3 * funq3 + funq1 / (kappa - 1);
}

/**********************************************************************/
/*                   MÉTODO DE EULER                                */
/**********************************************************************/
/*
   Función euler():
   Implementa el método de Euler para avanzar la solución en el tiempo.
   Se utiliza un esquema predictor-corrector:
     - Paso predictor (Euler simple) para obtener una solución provisional.
     - Paso corrector que promedia las cantidades conservadas para refinar la solución.
   
   Parámetros:
   - u1, u2, u3: Vectores para la solución actual (n+1).
   - u1old, u2old, u3old: Vectores para la solución del paso anterior (n).
   - u1tm1, u2tm1, u3tm1: Vectores para la solución del paso n-1 (usados para viscosidad artificial).
   - f1, f2, f3: Vectores para almacenar los flujos actuales.
   - f1pv, f2pv, f3pv: Vectores para almacenar los flujos del paso predictor.
   - q1v, q2v, q3v: Vectores para almacenar las cantidades conservadas actuales.
   - q1tv, q2tv, q3tv: Vectores para almacenar las cantidades conservadas del paso predictor.
   - u1tv, u2tv, u3tv: Vectores para almacenar las variables primitivas del paso predictor.
*/
void euler(gsl_vector *u1, gsl_vector *u2, gsl_vector *u3, 
           gsl_vector *u1old, gsl_vector *u2old, gsl_vector *u3old, 
           gsl_vector *u1tm1, gsl_vector *u2tm1, gsl_vector *u3tm1, 
           gsl_vector *f1, gsl_vector *f2, gsl_vector *f3, 
           gsl_vector *f1pv, gsl_vector *f2pv, gsl_vector *f3pv, 
           gsl_vector *q1v, gsl_vector *q2v, gsl_vector *q3v, 
           gsl_vector *q1tv, gsl_vector *q2tv, gsl_vector *q3tv, 
           gsl_vector *u1tv, gsl_vector *u2tv, gsl_vector *u3tv ) {

    double q1tempo, q2tempo, q3tempo, q1c, q2c, q3c;
    double u1c, u2c, u3c;
    double av1, av2, av3;
    
    for (int n = 1; n < NT; n++) {
        t = n * DT;
        openfile(t);

        /*************************************************************/
        /*                 CÁLCULO DE LOS FLUJOS                     */
        /*************************************************************/
        for (int i = 0; i < NX - 1; i++) {
            // Se calculan los flujos en cada celda usando la solución anterior
            fluxes(gsl_vector_get(u1old, i), 
                   gsl_vector_get(u2old, i), 
                   gsl_vector_get(u3old, i));
            gsl_vector_set(f1, i, funf1);
            gsl_vector_set(f2, i, funf2);
            gsl_vector_set(f3, i, funf3);
        }
        
        /*************************************************************/
        /*             PASO PREDICTOR (Euler Simple)                 */
        /*************************************************************/
        for (int i = 0; i < NX - 1; i++) {
            x = x_init + DX * i;
            if (i == 0) {
                // Para la frontera x = 0, se reestablecen las condiciones
                initialize(u1tv, u2tv, u3tv, t);
            } else {
                // Se calculan las cantidades conservadas a partir de la solución anterior
                charges(gsl_vector_get(u1old, i), 
                        gsl_vector_get(u2old, i), 
                        gsl_vector_get(u3old, i));
                
                // Se aplica el método de Euler para actualizar las cantidades
                q1tempo = q1 - (DT / DX) * (gsl_vector_get(f1, i) - gsl_vector_get(f1, i-1))
                          - (2 * DT / x) * gsl_vector_get(f1, i);
                q2tempo = q2 - (DT / DX) * (gsl_vector_get(f2, i) - gsl_vector_get(f2, i-1))
                          - (2 * DT / x) * gsl_vector_get(f2, i);
                q3tempo = q3 - (DT / DX) * (gsl_vector_get(f3, i) - gsl_vector_get(f3, i-1))
                          - (2 * DT / x) * gsl_vector_get(f3, i);
           
                // Se recuperan las variables primitivas a partir de las cantidades conservadas
                u1tempo = (q3tempo - 0.5 * q2tempo * q2tempo / q1tempo) * (kappa - 1);
                u2tempo = q1tempo;
                u3tempo = q2tempo / q1tempo;
           
                // Se almacenan los valores calculados en los vectores temporales
                gsl_vector_set(q1v, i, q1);
                gsl_vector_set(q2v, i, q2);
                gsl_vector_set(q3v, i, q3);
               
                gsl_vector_set(q1tv, i, q1tempo);
                gsl_vector_set(q2tv, i, q2tempo);
                gsl_vector_set(q3tv, i, q3tempo);
               
                gsl_vector_set(u1tv, i, u1tempo);
                gsl_vector_set(u2tv, i, u2tempo);
                gsl_vector_set(u3tv, i, u3tempo);
            }
        }
        
        /*************************************************************/
        /*         CÁLCULO DE LOS FLUJOS (PASO PREDICTOR)            */
        /*************************************************************/
        for (int i = 0; i < NX - 1; i++) {
            fluxes(gsl_vector_get(u1tv, i), 
                   gsl_vector_get(u2tv, i), 
                   gsl_vector_get(u3tv, i));
            gsl_vector_set(f1pv, i, funf1);
            gsl_vector_set(f2pv, i, funf2);
            gsl_vector_set(f3pv, i, funf3);
        }
        
        /*************************************************************/
        /*                   PASO CORRECTOR                        */
        /*************************************************************/
        for (int i = 0; i < NX - 1; i++) {
            x = x_init + DX * i;
            if (i == 0) {
                // Para la frontera x = 0 se reestablece la condición
                initialize(u1, u2, u3, t);
            } else {
                // Se promedian las cantidades conservadas y se corrige el paso
                q1c = 0.5 * ( gsl_vector_get(q1v, i) + gsl_vector_get(q1tv, i)
                         - (DT / DX) * (gsl_vector_get(f1pv, i) - gsl_vector_get(f1pv, i-1))
                         - (2 * DT / x) * gsl_vector_get(f1pv, i) );
                q2c = 0.5 * ( gsl_vector_get(q2v, i) + gsl_vector_get(q2tv, i)
                         - (DT / DX) * (gsl_vector_get(f2pv, i) - gsl_vector_get(f2pv, i-1))
                         - (2 * DT / x) * gsl_vector_get(f2pv, i) );
                q3c = 0.5 * ( gsl_vector_get(q3v, i) + gsl_vector_get(q3tv, i)
                         - (DT / DX) * (gsl_vector_get(f3pv, i) - gsl_vector_get(f3pv, i-1))
                         - (2 * DT / x) * gsl_vector_get(f3pv, i) );
           
                u1c = (q3c - 0.5 * q2c * q2c / q1c) * (kappa - 1);
                u2c = q1c;
                u3c = q2c / q1c;
              
                // Se aplica una viscosidad artificial para suavizar oscilaciones
                if (t != NT) {
                    av1 = (u1c - gsl_vector_get(u1old, i)) * (gsl_vector_get(u1old, i) - gsl_vector_get(u1tm1, i));
                    av2 = (u2c - gsl_vector_get(u2old, i)) * (gsl_vector_get(u2old, i) - gsl_vector_get(u2tm1, i));
                    av3 = (u3c - gsl_vector_get(u3old, i)) * (gsl_vector_get(u3old, i) - gsl_vector_get(u3tm1, i));
                  
                    if (av1 < 0) {
                        u1c = gsl_vector_get(u1old, i) + (1.0/4.0) * (u1c + gsl_vector_get(u1tm1, i) - 2 * gsl_vector_get(u1old, i));
                    }
                    if (av2 < 0) {
                        u2c = gsl_vector_get(u2old, i) + (1.0/4.0) * (u2c + gsl_vector_get(u2tm1, i) - 2 * gsl_vector_get(u2old, i));
                    }
                    if (av3 < 0) {
                        u3c = gsl_vector_get(u3old, i) + (1.0/4.0) * (u3c + gsl_vector_get(u3tm1, i) - 2 * gsl_vector_get(u3old, i));
                    }
                }
           
                gsl_vector_set(u1, i, u1c);
                gsl_vector_set(u2, i, u2c);  
                gsl_vector_set(u3, i, u3c);
            }
           
            // Se escribe la solución del paso actual en el archivo
            fprintf(file, "%g\t%g\t%g\t%g\n", x, 
                    gsl_vector_get(u1, i), 
                    gsl_vector_get(u2, i),
                    gsl_vector_get(u3, i));
        }

        // Actualizar vectores para el siguiente paso:
        // u_tm1 <- u_old y u_old <- u
        gsl_vector_memcpy(u1tm1, u1old);
        gsl_vector_memcpy(u2tm1, u2old);  
        gsl_vector_memcpy(u3tm1, u3old);
        
        gsl_vector_memcpy(u1old, u1);
        gsl_vector_memcpy(u2old, u2);  
        gsl_vector_memcpy(u3old, u3);
        
        closefile();
    }
}

/**********************************************************************/
/*                              MAIN                                  */
/**********************************************************************/
/*
   Función main():
   - Inicializa los parámetros y reserva memoria para los vectores.
   - Establece las condiciones iniciales (suavizadas) mediante initialize().
   - Ejecuta el método de Euler para evolucionar la solución en el tiempo.
   - Libera la memoria reservada y termina la simulación.
*/
int main() {
    // Inicialización de parámetros globales
    initialize_variables();
    
    // Reserva de memoria para las variables primitivas:
    // u1: presión, u2: densidad, u3: velocidad.
    gsl_vector *u1 = gsl_vector_alloc(NX + 1);   
    gsl_vector *u2 = gsl_vector_alloc(NX + 1); 
    gsl_vector *u3 = gsl_vector_alloc(NX + 1); 
    
    // Vectores para la solución del paso anterior (n)
    gsl_vector *u1old = gsl_vector_alloc(NX + 1); 
    gsl_vector *u2old = gsl_vector_alloc(NX + 1);
    gsl_vector *u3old = gsl_vector_alloc(NX + 1);
    
    // Vectores para la solución en el paso n-1 (para viscosidad artificial)
    gsl_vector *u1tm1 = gsl_vector_alloc(NX + 1); 
    gsl_vector *u2tm1 = gsl_vector_alloc(NX + 1);
    gsl_vector *u3tm1 = gsl_vector_alloc(NX + 1);
    
    // Vectores para almacenar los flujos actuales
    gsl_vector *f1 = gsl_vector_alloc(NX + 1);   
    gsl_vector *f2 = gsl_vector_alloc(NX + 1);
    gsl_vector *f3 = gsl_vector_alloc(NX + 1);
    
    // Vectores para almacenar los flujos del paso predictor
    gsl_vector *f1pv = gsl_vector_alloc(NX + 1);   
    gsl_vector *f2pv = gsl_vector_alloc(NX + 1); 
    gsl_vector *f3pv = gsl_vector_alloc(NX + 1);
    
    // Vectores para almacenar las cantidades conservadas actuales
    gsl_vector *q1v = gsl_vector_alloc(NX + 1);   
    gsl_vector *q2v = gsl_vector_alloc(NX + 1); 
    gsl_vector *q3v = gsl_vector_alloc(NX + 1);
    
    // Vectores para almacenar las cantidades conservadas del paso predictor
    gsl_vector *q1tv = gsl_vector_alloc(NX + 1);   
    gsl_vector *q2tv = gsl_vector_alloc(NX + 1); 
    gsl_vector *q3tv = gsl_vector_alloc(NX + 1);
    
    // Vectores para almacenar las variables primitivas del paso predictor
    gsl_vector *u1tv = gsl_vector_alloc(NX + 1);   
    gsl_vector *u2tv = gsl_vector_alloc(NX + 1); 
    gsl_vector *u3tv = gsl_vector_alloc(NX + 1);
    
    // Inicializa la malla y establece las condiciones iniciales suavizadas
    initialize(u1old, u2old, u3old, t_init);
    
    // Copia inicial para la solución en el paso n-1
    gsl_vector_memcpy(u1tm1, u1old);
    gsl_vector_memcpy(u2tm1, u2old);
    gsl_vector_memcpy(u3tm1, u3old);

    // Ejecuta el método de Euler para la evolución temporal
    euler(u1, u2, u3, 
          u1old, u2old, u3old, 
          u1tm1, u2tm1, u3tm1, 
          f1, f2, f3, 
          f1pv, f2pv, f3pv, 
          q1v, q2v, q3v, 
          q1tv, q2tv, q3tv, 
          u1tv, u2tv, u3tv );

    // Liberación de la memoria reservada para todos los vectores
    gsl_vector_free(u2);
    gsl_vector_free(u3);
    gsl_vector_free(u1);
    
    gsl_vector_free(u2old);
    gsl_vector_free(u3old);
    gsl_vector_free(u1old);
    
    gsl_vector_free(u2tm1);
    gsl_vector_free(u3tm1);
    gsl_vector_free(u1tm1);
    
    gsl_vector_free(f2);
    gsl_vector_free(f3);
    gsl_vector_free(f1);
    
    gsl_vector_free(f1pv);
    gsl_vector_free(f2pv);
    gsl_vector_free(f3pv);
    
    gsl_vector_free(q1v);
    gsl_vector_free(q2v);
    gsl_vector_free(q3v);
    
    gsl_vector_free(q1tv);
    gsl_vector_free(q2tv);
    gsl_vector_free(q3tv);

    return 0;
}

/**********************************************************************/
/*                  FUNCIONES DE MANEJO DE ARCHIVOS                   */
/**********************************************************************/
/*
   Función openfile():
   - Establece el nombre del archivo utilizando el valor del tiempo.
   - Llama a writeToFile() para abrir el archivo en modo escritura.
*/
void openfile(double floatValue) {
    snprintf(filename, sizeof(filename), "data_%.5f.dat", floatValue);
    writeToFile(floatValue);
}

/*
   Función closefile():
   - Cierra el archivo actualmente abierto.
*/
void closefile() {
    fclose(file);
}

/*
   Función writeToFile():
   - Abre el archivo cuyo nombre está en la variable filename en modo escritura.
   - Si falla la apertura, muestra un mensaje de error.
*/
void writeToFile(double value) {
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error abriendo el archivo");
        return;
    }
}

