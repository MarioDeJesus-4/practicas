#include <stdio.h>
#include <stdlib.h>
#include <math.h>           // Se requiere para tanh()
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

/****************************************************/
/*  This code solves a 1D fluid flow problem        */
/*  (presumably the Euler equations) using a        */
/*  finite-difference approach and a 4th-order      */
/*  Runge-Kutta (RK4) time integration scheme.      */
/****************************************************/

/* 
   File management variables and functions.
*/
void writeToFile(double);
void openfile(double);
void closefile();
char filename[50];  // File name
FILE *file;         // File pointer

/* 
   Global simulation parameters.
*/
double x_final, x_init, NX, DX;
double t_final, t_init, NT, DT;
double kappa;

/*
   Auxiliary variables.
*/
double t, x;
double u2tempo, u3tempo, u1tempo;
double u1i, u2i, u3i;
double funf1, funf2, funf3;
double q1, q2, q3;

/*
   Initializes the global simulation parameters.
*/
void initialize_variables() {
    x_final = 10.0;
    x_init  = 1.0;
    NX = 1000; 
    DX = (x_final - x_init) / NX; // Spatial step size

    t_final = 1.0;
    t_init  = 0.0;
    NT = 1000; 
    DT = (t_final - t_init) / NT; // Time step size

    kappa = 1.4;  // Adiabatic coefficient
}

/*
   Initializes the grid with the initial conditions.
   The initial conditions are now defined with a smooth transition 
   using a tanh function, instead of a sharp discontinuity.
   
   - u1in, u2in, u3in: vectors to store pressure, density, and velocity.
   - tin: current time. If tin==0, initial conditions are set.
   
   The left state is defined as:
       u1_left = 5.0, u2_left = 4.0, u3_left = 5.75,
   and the right state as:
       u1_right = 0.5, u2_right = 0.5, u3_right = 0.75.
   The smooth transition is centered at x0 with width controlled by epsilon.
*/
void initialize(gsl_vector *u1in, gsl_vector *u2in, gsl_vector *u3in, double tin) {
    if (tin == 0) {
        openfile(tin); 

        // Define center of the domain and epsilon for the smooth transition
        double x0 = (x_init + x_final) / 2.0;
        double epsilon = (x_final - x_init) / 50.0; // Ajustar este parámetro según se requiera

        for (int i = 0; i < NX; i++) {
            x = x_init + DX * i;
            
            // Smooth transition using tanh:
            // u = u_left + 0.5*(u_right - u_left)*(1 + tanh((x - x0)/epsilon))
            u1tempo = 5.0 + 0.5 * (0.5 - 5.0) * (1 + tanh((x - x0) / epsilon));
            u2tempo = 4.0 + 0.5 * (0.5 - 4.0) * (1 + tanh((x - x0) / epsilon));
            u3tempo = 5.75 + 0.5 * (0.75 - 5.75) * (1 + tanh((x - x0) / epsilon));

            gsl_vector_set(u1in, i, u1tempo);
            gsl_vector_set(u2in, i, u2tempo);
            gsl_vector_set(u3in, i, u3tempo);
            
            // Write initial condition for the current cell
            fprintf(file, "%g\t%g\t%g\t%g\n", x,
                    gsl_vector_get(u1in, i),
                    gsl_vector_get(u2in, i),
                    gsl_vector_get(u3in, i));
            
            // Save condition at x = 0 for boundary purposes.
            if (i == 0) {
                u1i = gsl_vector_get(u1in, 0);
                u2i = gsl_vector_get(u2in, 0);
                u3i = gsl_vector_get(u3in, 0);
            }
        }
        closefile();
    } else {
        // For subsequent times, enforce the boundary condition at x=0.
        gsl_vector_set(u1in, 0, u1i);  
        gsl_vector_set(u2in, 0, u2i);
        gsl_vector_set(u3in, 0, u3i);
    }
}

/*
   Computes the fluxes (f1, f2, f3) given u1, u2, and u3.
   Assumes: u1 = pressure, u2 = density, u3 = velocity.
*/
void fluxes(double funu1, double funu2, double funu3) {
    funf1 = funu2 * funu3;
    funf2 = funu2 * funu3 * funu3 + funu1;
    funf3 = funu2 * funu3 * (
        0.5 * funu3 * funu3 + (kappa / (kappa - 1)) * funu1 / funu2
    );
}

/*
   Computes the conserved quantities (q1, q2, q3) from (u1, u2, u3).
   For instance, q1 = density, q2 = momentum, q3 = total energy.
*/
void charges(double funq1, double funq2, double funq3) {
    q1 = funq2;              // density
    q2 = funq2 * funq3;      // momentum
    q3 = 0.5 * funq2 * funq3 * funq3 + funq1 / (kappa - 1);
}

/*
   Computes the time derivative of Q (dQ/dt) based on the flux differences 
   and a geometric term (2*F/x). The derivative is stored in dq1, dq2, dq3.
*/
void compute_time_derivative(gsl_vector *u1old, gsl_vector *u2old, gsl_vector *u3old,
                             gsl_vector *dq1, gsl_vector *dq2, gsl_vector *dq3,
                             int n, int totalNT, double time) 
{
    gsl_vector *f1_local = gsl_vector_alloc(NX + 1);
    gsl_vector *f2_local = gsl_vector_alloc(NX + 1);
    gsl_vector *f3_local = gsl_vector_alloc(NX + 1);

    // Refresh boundary condition at x=0
    initialize(u1old, u2old, u3old, time);

    // Compute fluxes in each cell
    for (int i = 0; i < NX; i++) {
        fluxes(gsl_vector_get(u1old, i),
               gsl_vector_get(u2old, i),
               gsl_vector_get(u3old, i));
        gsl_vector_set(f1_local, i, funf1);
        gsl_vector_set(f2_local, i, funf2);
        gsl_vector_set(f3_local, i, funf3);
    }

    // Compute dQ/dt = -[F(i)-F(i-1)]/DX - (2*F(i))/x, for i >= 1
    for (int i = 0; i < NX; i++) {
        if (i == 0) {
            gsl_vector_set(dq1, i, 0.0);
            gsl_vector_set(dq2, i, 0.0);
            gsl_vector_set(dq3, i, 0.0);
        } else {
            x = x_init + DX * i;
            charges(gsl_vector_get(u1old, i),
                    gsl_vector_get(u2old, i),
                    gsl_vector_get(u3old, i));
            
            double flux_im1_f1 = gsl_vector_get(f1_local, i) - gsl_vector_get(f1_local, i - 1);
            double flux_im1_f2 = gsl_vector_get(f2_local, i) - gsl_vector_get(f2_local, i - 1);
            double flux_im1_f3 = gsl_vector_get(f3_local, i) - gsl_vector_get(f3_local, i - 1);

            double dq1_dt = - (flux_im1_f1 / DX) - (2.0 / x) * gsl_vector_get(f1_local, i);
            double dq2_dt = - (flux_im1_f2 / DX) - (2.0 / x) * gsl_vector_get(f2_local, i);
            double dq3_dt = - (flux_im1_f3 / DX) - (2.0 / x) * gsl_vector_get(f3_local, i);

            gsl_vector_set(dq1, i, dq1_dt);
            gsl_vector_set(dq2, i, dq2_dt);
            gsl_vector_set(dq3, i, dq3_dt);
        }
    }
    
    gsl_vector_free(f1_local);
    gsl_vector_free(f2_local);
    gsl_vector_free(f3_local);
}

/*
   Implements the 4th-order Runge-Kutta (RK4) time integration method.
   - u1, u2, u3: vectors with the solution at the new time step (n+1)
   - u1old, u2old, u3old: vectors with the solution at the previous time step (n)
   - u1tm1, u2tm1, u3tm1: vectors with the solution at time step (n-1) used in artificial viscosity.
*/
void runge_kutta_4(gsl_vector *u1, gsl_vector *u2, gsl_vector *u3,
                   gsl_vector *u1old, gsl_vector *u2old, gsl_vector *u3old,
                   gsl_vector *u1tm1, gsl_vector *u2tm1, gsl_vector *u3tm1)
{
    // --- 1) Allocate temporary vectors ---
    gsl_vector *k1_q1 = gsl_vector_alloc(NX + 1);
    gsl_vector *k1_q2 = gsl_vector_alloc(NX + 1);
    gsl_vector *k1_q3 = gsl_vector_alloc(NX + 1);

    gsl_vector *k2_q1 = gsl_vector_alloc(NX + 1);
    gsl_vector *k2_q2 = gsl_vector_alloc(NX + 1);
    gsl_vector *k2_q3 = gsl_vector_alloc(NX + 1);

    gsl_vector *k3_q1 = gsl_vector_alloc(NX + 1);
    gsl_vector *k3_q2 = gsl_vector_alloc(NX + 1);
    gsl_vector *k3_q3 = gsl_vector_alloc(NX + 1);

    gsl_vector *k4_q1 = gsl_vector_alloc(NX + 1);
    gsl_vector *k4_q2 = gsl_vector_alloc(NX + 1);
    gsl_vector *k4_q3 = gsl_vector_alloc(NX + 1);

    gsl_vector *q1_temp = gsl_vector_alloc(NX + 1);
    gsl_vector *q2_temp = gsl_vector_alloc(NX + 1);
    gsl_vector *q3_temp = gsl_vector_alloc(NX + 1);

    // Intermediate vectors for the half step solution
    gsl_vector *u1_half = gsl_vector_alloc(NX + 1);
    gsl_vector *u2_half = gsl_vector_alloc(NX + 1);
    gsl_vector *u3_half = gsl_vector_alloc(NX + 1);

    // --- 2) Time loop ---
    for (int n = 1; n <= NT; n++) {
        t = n * DT;

        // -------- k1 ----------
        compute_time_derivative(u1old, u2old, u3old,
                                k1_q1, k1_q2, k1_q3,
                                n, (int)NT, t - DT);

        // Q^n + (DT/2)*k1
        for (int i = 0; i < NX; i++) {
            charges(gsl_vector_get(u1old, i),
                    gsl_vector_get(u2old, i),
                    gsl_vector_get(u3old, i));
            double Q1n = q1;
            double Q2n = q2;
            double Q3n = q3;

            double Q1n_half = Q1n + 0.5 * DT * gsl_vector_get(k1_q1, i);
            double Q2n_half = Q2n + 0.5 * DT * gsl_vector_get(k1_q2, i);
            double Q3n_half = Q3n + 0.5 * DT * gsl_vector_get(k1_q3, i);

            gsl_vector_set(q1_temp, i, Q1n_half);
            gsl_vector_set(q2_temp, i, Q2n_half);
            gsl_vector_set(q3_temp, i, Q3n_half);
        }

        // Convert Q^(n+1/2) to u^(n+1/2)
        for (int i = 0; i < NX; i++) {
            double Q1_val = gsl_vector_get(q1_temp, i);
            double Q2_val = gsl_vector_get(q2_temp, i);
            double Q3_val = gsl_vector_get(q3_temp, i);

            double new_u1 = (Q3_val - 0.5 * (Q2_val * Q2_val) / Q1_val) * (kappa - 1.0);
            double new_u2 = Q1_val;
            double new_u3 = Q2_val / Q1_val;

            gsl_vector_set(u1_half, i, new_u1);
            gsl_vector_set(u2_half, i, new_u2);
            gsl_vector_set(u3_half, i, new_u3);
        }

        // -------- k2 ----------
        compute_time_derivative(u1_half, u2_half, u3_half,
                                k2_q1, k2_q2, k2_q3,
                                n, (int)NT, t - (DT / 2.0));

        for (int i = 0; i < NX; i++) {
            charges(gsl_vector_get(u1old, i),
                    gsl_vector_get(u2old, i),
                    gsl_vector_get(u3old, i));
            double Q1n = q1;
            double Q2n = q2;
            double Q3n = q3;

            double Q1n_half = Q1n + 0.5 * DT * gsl_vector_get(k2_q1, i);
            double Q2n_half = Q2n + 0.5 * DT * gsl_vector_get(k2_q2, i);
            double Q3n_half = Q3n + 0.5 * DT * gsl_vector_get(k2_q3, i);

            gsl_vector_set(q1_temp, i, Q1n_half);
            gsl_vector_set(q2_temp, i, Q2n_half);
            gsl_vector_set(q3_temp, i, Q3n_half);
        }

        for (int i = 0; i < NX; i++) {
            double Q1_val = gsl_vector_get(q1_temp, i);
            double Q2_val = gsl_vector_get(q2_temp, i);
            double Q3_val = gsl_vector_get(q3_temp, i);

            double new_u1 = (Q3_val - 0.5 * (Q2_val * Q2_val) / Q1_val) * (kappa - 1.0);
            double new_u2 = Q1_val;
            double new_u3 = Q2_val / Q1_val;

            gsl_vector_set(u1_half, i, new_u1);
            gsl_vector_set(u2_half, i, new_u2);
            gsl_vector_set(u3_half, i, new_u3);
        }

        // -------- k3 ----------
        compute_time_derivative(u1_half, u2_half, u3_half,
                                k3_q1, k3_q2, k3_q3,
                                n, (int)NT, t - (DT / 2.0));

        for (int i = 0; i < NX; i++) {
            charges(gsl_vector_get(u1old, i),
                    gsl_vector_get(u2old, i),
                    gsl_vector_get(u3old, i));
            double Q1n = q1;
            double Q2n = q2;
            double Q3n = q3;

            double Q1n_1 = Q1n + DT * gsl_vector_get(k3_q1, i);
            double Q2n_1 = Q2n + DT * gsl_vector_get(k3_q2, i);
            double Q3n_1 = Q3n + DT * gsl_vector_get(k3_q3, i);

            gsl_vector_set(q1_temp, i, Q1n_1);
            gsl_vector_set(q2_temp, i, Q2n_1);
            gsl_vector_set(q3_temp, i, Q3n_1);
        }

        for (int i = 0; i < NX; i++) {
            double Q1_val = gsl_vector_get(q1_temp, i);
            double Q2_val = gsl_vector_get(q2_temp, i);
            double Q3_val = gsl_vector_get(q3_temp, i);

            double new_u1 = (Q3_val - 0.5 * (Q2_val * Q2_val) / Q1_val) * (kappa - 1.0);
            double new_u2 = Q1_val;
            double new_u3 = Q2_val / Q1_val;

            gsl_vector_set(u1_half, i, new_u1);
            gsl_vector_set(u2_half, i, new_u2);
            gsl_vector_set(u3_half, i, new_u3);
        }

        // -------- k4 ----------
        compute_time_derivative(u1_half, u2_half, u3_half,
                                k4_q1, k4_q2, k4_q3,
                                n, (int)NT, t);

        // -------- Combine k1, k2, k3, k4 to obtain Q^(n+1) and then u^(n+1) --------
        for (int i = 0; i < NX; i++) {
            charges(gsl_vector_get(u1old, i),
                    gsl_vector_get(u2old, i),
                    gsl_vector_get(u3old, i));
            double Q1n = q1;
            double Q2n = q2;
            double Q3n = q3;

            double sum_k1 = gsl_vector_get(k1_q1, i);
            double sum_k2 = gsl_vector_get(k2_q1, i);
            double sum_k3 = gsl_vector_get(k3_q1, i);
            double sum_k4 = gsl_vector_get(k4_q1, i);
            double Q1n_next = Q1n + (DT / 6.0) * (sum_k1 + 2.0 * sum_k2 + 2.0 * sum_k3 + sum_k4);

            sum_k1 = gsl_vector_get(k1_q2, i);
            sum_k2 = gsl_vector_get(k2_q2, i);
            sum_k3 = gsl_vector_get(k3_q2, i);
            sum_k4 = gsl_vector_get(k4_q2, i);
            double Q2n_next = Q2n + (DT / 6.0) * (sum_k1 + 2.0 * sum_k2 + 2.0 * sum_k3 + sum_k4);

            sum_k1 = gsl_vector_get(k1_q3, i);
            sum_k2 = gsl_vector_get(k2_q3, i);
            sum_k3 = gsl_vector_get(k3_q3, i);
            sum_k4 = gsl_vector_get(k4_q3, i);
            double Q3n_next = Q3n + (DT / 6.0) * (sum_k1 + 2.0 * sum_k2 + 2.0 * sum_k3 + sum_k4);

            // Convert Q^(n+1) -> u^(n+1)
            double new_u1 = (Q3n_next - 0.5 * (Q2n_next * Q2n_next) / Q1n_next) * (kappa - 1.0);
            double new_u2 = Q1n_next;
            double new_u3 = Q2n_next / Q1n_next;

            // --- Artificial Viscosity ---
            if (n != (int)NT) {
                double av1 = (new_u1 - gsl_vector_get(u1old, i)) *
                             (gsl_vector_get(u1old, i) - gsl_vector_get(u1tm1, i));
                double av2 = (new_u2 - gsl_vector_get(u2old, i)) *
                             (gsl_vector_get(u2old, i) - gsl_vector_get(u2tm1, i));
                double av3 = (new_u3 - gsl_vector_get(u3old, i)) *
                             (gsl_vector_get(u3old, i) - gsl_vector_get(u3tm1, i));
                
                if (av1 < 0.0) {
                    new_u1 = gsl_vector_get(u1old, i) +
                             0.25 * (new_u1 + gsl_vector_get(u1tm1, i) - 2.0 * gsl_vector_get(u1old, i));
                }
                if (av2 < 0.0) {
                    new_u2 = gsl_vector_get(u2old, i) +
                             0.25 * (new_u2 + gsl_vector_get(u2tm1, i) - 2.0 * gsl_vector_get(u2old, i));
                }
                if (av3 < 0.0) {
                    new_u3 = gsl_vector_get(u3old, i) +
                             0.25 * (new_u3 + gsl_vector_get(u3tm1, i) - 2.0 * gsl_vector_get(u3old, i));
                }
            }

            gsl_vector_set(u1, i, new_u1);
            gsl_vector_set(u2, i, new_u2);
            gsl_vector_set(u3, i, new_u3);
        }

        // Write the solution for time t to file
        openfile(t);
        for (int i = 0; i < NX; i++) {
            x = x_init + DX * i;
            fprintf(file, "%g\t%g\t%g\t%g\n", x,
                    gsl_vector_get(u1, i),
                    gsl_vector_get(u2, i),
                    gsl_vector_get(u3, i));
        }
        closefile();

        // Update: u_tm1 <- u_old and u_old <- u
        gsl_vector_memcpy(u1tm1, u1old);
        gsl_vector_memcpy(u2tm1, u2old);
        gsl_vector_memcpy(u3tm1, u3old);

        gsl_vector_memcpy(u1old, u1);
        gsl_vector_memcpy(u2old, u2);
        gsl_vector_memcpy(u3old, u3);
    }

    // --- 3) Free temporary vectors ---
    gsl_vector_free(k1_q1);
    gsl_vector_free(k1_q2);
    gsl_vector_free(k1_q3);

    gsl_vector_free(k2_q1);
    gsl_vector_free(k2_q2);
    gsl_vector_free(k2_q3);

    gsl_vector_free(k3_q1);
    gsl_vector_free(k3_q2);
    gsl_vector_free(k3_q3);

    gsl_vector_free(k4_q1);
    gsl_vector_free(k4_q2);
    gsl_vector_free(k4_q3);

    gsl_vector_free(q1_temp);
    gsl_vector_free(q2_temp);
    gsl_vector_free(q3_temp);

    gsl_vector_free(u1_half);
    gsl_vector_free(u2_half);
    gsl_vector_free(u3_half);
}

/*
   Main function.
*/
int main() {
    // Initialize global parameters
    initialize_variables();
    
    // Allocate memory for vectors
    gsl_vector *u1 = gsl_vector_alloc(NX + 1);   
    gsl_vector *u2 = gsl_vector_alloc(NX + 1); 
    gsl_vector *u3 = gsl_vector_alloc(NX + 1); 
    
    gsl_vector *u1old = gsl_vector_alloc(NX + 1); 
    gsl_vector *u2old = gsl_vector_alloc(NX + 1);
    gsl_vector *u3old = gsl_vector_alloc(NX + 1);
    
    gsl_vector *u1tm1 = gsl_vector_alloc(NX + 1); 
    gsl_vector *u2tm1 = gsl_vector_alloc(NX + 1);
    gsl_vector *u3tm1 = gsl_vector_alloc(NX + 1);
    
    // Initialize grid and initial conditions
    initialize(u1old, u2old, u3old, t_init);
    
    // Copy initial conditions for step n-1
    gsl_vector_memcpy(u1tm1, u1old);
    gsl_vector_memcpy(u2tm1, u2old);
    gsl_vector_memcpy(u3tm1, u3old);

    // Run RK4 method
    runge_kutta_4(u1, u2, u3,
                  u1old, u2old, u3old,
                  u1tm1, u2tm1, u3tm1);

    // Free allocated memory
    gsl_vector_free(u2);
    gsl_vector_free(u3);
    gsl_vector_free(u1);
    
    gsl_vector_free(u2old);
    gsl_vector_free(u3old);
    gsl_vector_free(u1old);
    
    gsl_vector_free(u2tm1);
    gsl_vector_free(u3tm1);
    gsl_vector_free(u1tm1);

    return 0;
}

/*
   File handling functions.
*/
void openfile(double floatValue) {
    snprintf(filename, sizeof(filename), "data_%.5f.dat", floatValue);
    writeToFile(floatValue);  
}

void closefile() {
    fclose(file); 
}

void writeToFile(double value) {
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
}

