#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;
double dt, dt_old; /* Alterado de static para global */

double Random(void)
    /* ----------------------------------------------------------------
     * Random returns a pseudo-random real number uniformly distributed
     * between 0.0 and 1.0.
     * ----------------------------------------------------------------
     */
{
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;
    long t;

    t = MULTIPLIER * (seed % Q) - R * (seed / Q);

    if (t > 0) {
        seed = t;
    } else {
        seed = t + MODULUS;
    }

    return (double) seed / MODULUS;
}

/*
 * End of the pRNG algorithm
 */

typedef struct {
    double x, y, z;
    double mass;
} Particle;

typedef struct {
    double xold, yold, zold;
    double fx, fy, fz;
} ParticleV;

void InitParticles(Particle[], ParticleV [], int , int);
double ComputeForces(Particle [], Particle [], ParticleV [], int, int, int);
double ComputeNewPos(Particle [], ParticleV [], int, int, double);

int start(int rank, double ppp)
{
    return (int)(ppp*(rank-1));
}

int end(int rank, double ppp)
{
    return (int)(ppp*(rank));
}

typedef enum {
    COMPUTE_FORCES,
    COMPUTE_POS,
    STOP,
} Instruction;

typedef enum {
    INSTRUCTION,
    MAX_F,
    PARTICLES,
    PV,
    RANK,
} Tags;

void master(int npart, int nsteps, int nproc)
{
    int i, s, e;
    int src;
    Instruction inst;
    MPI_Status st;
    double ppp;
    double max_f, local_max_f;

    Particle *particles;
    ParticleV *pv;

    particles = malloc(sizeof(Particle)*npart);
    pv = malloc(sizeof(ParticleV)*npart);
    ppp = npart/nproc;

    InitParticles(particles, pv, 0, npart);

    while (--nsteps) {
        inst = COMPUTE_FORCES;

        for (i = 0; i < npart; i++) {
            fprintf(stdout, "%.5lf %.5lf\n", particles[i].x, particles[i].y);
        }

        for (i = 1; i < nproc; i++) {
            MPI_Send(&inst, 1, MPI_INT, i, INSTRUCTION, MPI_COMM_WORLD);
            MPI_Send(particles, sizeof(Particle)*npart, MPI_BYTE, i, PARTICLES,
                    MPI_COMM_WORLD);
            MPI_Send(pv, sizeof(ParticleV)*npart, MPI_BYTE, i, PV,
                    MPI_COMM_WORLD);
        }

        for (i = 1; i < nproc; i++) {
            MPI_Recv(&local_max_f, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MAX_F,
                    MPI_COMM_WORLD, &st);
            max_f = local_max_f > max_f ? local_max_f : max_f;
            src = st.MPI_SOURCE;
            s = start(src, ppp);
            e = end(src, ppp);
            MPI_Recv(&pv[s], sizeof(ParticleV)*(e-s), MPI_BYTE, src,
                    PV, MPI_COMM_WORLD, &st);
        }

        inst = COMPUTE_POS;

        for (i = 1; i < nproc; i++) {
            MPI_Send(&inst, 1, MPI_INT, i, INSTRUCTION, MPI_COMM_WORLD);
            MPI_Send(&max_f, 1, MPI_DOUBLE, i, MAX_F, MPI_COMM_WORLD);
            MPI_Send(pv, sizeof(ParticleV)*npart, MPI_BYTE, i, PV,
                    MPI_COMM_WORLD);
        }

        for (i = 1; i < nproc; i++) {
            MPI_Recv(&src, 1, MPI_INT, MPI_ANY_SOURCE, RANK, MPI_COMM_WORLD, &st);
            s = start(src, ppp);
            e = end(src, ppp);
            MPI_Recv(&pv[s], sizeof(ParticleV)*(e-s), MPI_BYTE, src,
                    PV, MPI_COMM_WORLD, &st);
            MPI_Recv(&particles[s], sizeof(Particle)*(e-s), MPI_BYTE, src,
                    PARTICLES, MPI_COMM_WORLD, &st);
        }
    }

    for (i = 1; i < nproc; i++) {
        inst = STOP;
        MPI_Send(&inst, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    for (i = 0; i < npart; i++) {
        fprintf(stdout, "%.5lf %.5lf\n", particles[i].x, particles[i].y);
    }
}

void slave(int rank, int npart, int nproc)
{
    Instruction instruction;
    MPI_Status st;
    int running;
    int i, s, e;
    Particle *particles;
    ParticleV *pv;
    double ppp;
    double max_f, sim_t;

    particles = malloc(sizeof(Particle)*npart);
    pv = malloc(sizeof(ParticleV)*npart);

    ppp = (double)npart/nproc;

    s = start(rank, ppp);
    e = end(rank, ppp);

    running = 1;

    while (running) {
        MPI_Recv(&instruction, 1, MPI_INT, 0, INSTRUCTION, MPI_COMM_WORLD, &st);

        switch (instruction) {
            case COMPUTE_FORCES:
                MPI_Recv(particles, sizeof(Particle)*npart, MPI_BYTE, 0,
                        PARTICLES, MPI_COMM_WORLD, &st);
                MPI_Recv(pv, sizeof(ParticleV)*npart, MPI_BYTE, 0,
                        PV, MPI_COMM_WORLD, &st);
                for (i = 0; i < npart; i++) {
                    fprintf(stdout, "%d: %.5lf %.5lf\n", rank, particles[i].x, particles[i].y);
                }
                max_f = ComputeForces(particles, particles, pv, s, e, npart);
                MPI_Send(&max_f, 1, MPI_DOUBLE, 0, MAX_F, MPI_COMM_WORLD);
                MPI_Send(&pv[s], sizeof(ParticleV)*(e-s), MPI_BYTE, 0, PV,
                        MPI_COMM_WORLD);
                break;
            case COMPUTE_POS:
                MPI_Recv(pv, sizeof(ParticleV)*npart, MPI_BYTE, 0,
                        PV, MPI_COMM_WORLD, &st);
                MPI_Recv(&max_f, 1, MPI_DOUBLE, 0,
                        MAX_F, MPI_COMM_WORLD, &st);
                sim_t += ComputeNewPos(particles, pv, s, e, max_f);
                MPI_Send(&rank, 1, MPI_INT, 0, RANK, MPI_COMM_WORLD);
                MPI_Send(&pv[s], sizeof(ParticleV)*(e-s), MPI_BYTE, 0, PV,
                        MPI_COMM_WORLD);
                MPI_Send(&particles[s], sizeof(Particle)*(e-s), MPI_BYTE, 0,
                        PARTICLES, MPI_COMM_WORLD);
                break;
            case STOP:
                running = 0;
        }
    }
}

int main(int argc, char **argv)
{
    int rank, nproc;
    int npart;
    int nsteps;                 /* number of times in loop */
    double sim_t;            /* Simulation time */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 3) {
        if (rank == 0) {
            printf("Wrong number of parameters.\n"
                    "Usage: nbody num_bodies timesteps\n");
        }
        exit(1);
    }

    npart = atoi(argv[1]);
    nsteps = atoi(argv[2]);

    if (rank == 0) {
        master(npart, nsteps, nproc);
    } else {
        slave(rank, nproc, npart);
    }

    MPI_Finalize();

    // sim_t = 0.0;

    return 0;
}

void InitParticles(Particle particles[], ParticleV pv[], int start, int end)
{
    int i;
    for (i = start; i < end; i++) {
        particles[i].x    = Random();
        particles[i].y    = Random();
        particles[i].z    = Random();
        particles[i].mass = 1.0;
        pv[i].xold    = particles[i].x;
        pv[i].yold    = particles[i].y;
        pv[i].zold    = particles[i].z;
        pv[i].fx      = 0;
        pv[i].fy      = 0;
        pv[i].fz      = 0;
    }
}

double ComputeForces(Particle myparticles[], Particle others[], ParticleV pv[], int start, int end, int npart)
{
    double max_f;
    int i;

    max_f = 0.0;

    for (i = start; i < end; i++) {
        int j;
        double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
        rmin = 100.0;
        xi   = myparticles[i].x;
        yi   = myparticles[i].y;
        fx   = 0.0;
        fy   = 0.0;

        for (j = 0; j < npart; j++) {
            rx = xi - others[j].x;
            ry = yi - others[j].y;
            mj = others[j].mass;
            r  = rx * rx + ry * ry;
            /* ignore overlap and same particle */
            if (r == 0.0) {continue; }
            if (r < rmin) {rmin = r; }
            r  = r * sqrt(r);
            fx -= mj * rx / r;
            fy -= mj * ry / r;
        }

        pv[i].fx += fx;
        pv[i].fy += fy;
        fx = sqrt(fx*fx + fy*fy)/rmin;
        if (fx > max_f) {max_f = fx; }
    }

    return max_f;
}

double ComputeNewPos(Particle particles[], ParticleV pv[], int start, int end, double max_f)
{
    int i;
    double a0, a1, a2;
    double dt_new;

    a0   = 2.0 / (dt * (dt + dt_old));
    a2   = 2.0 / (dt_old * (dt + dt_old));
    a1   = -(a0 + a2);

    for (i = start; i < end; i++) {
        double xi, yi;
        xi             = particles[i].x;
        yi             = particles[i].y;
        particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
        particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
        pv[i].xold     = xi;
        pv[i].yold     = yi;
        pv[i].fx       = 0;
        pv[i].fy       = 0;
    }

    dt_new = 1.0/sqrt(max_f);

    /* Set a minimum: */
    if (dt_new < 1.0e-6) {dt_new = 1.0e-6; }

    /* Modify time step */
    if (dt_new < dt) {
        dt_old = dt;
        dt     = dt_new;
    } else if (dt_new > 4.0 * dt) {
        dt_old = dt;
        dt    *= 2.0;
    }

    return dt_old;
}

