#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;
double dt, dt_old; /* Alterado de static para global */
pthread_mutex_t mutex;

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
  if (t > 0)
    seed = t;
  else
    seed = t + MODULUS;
  return ((double) seed / MODULUS);
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

typedef struct {
	Particle ps;
	ParticleV psv;
}	PnV;

void InitParticles(int );
//void *InitParticles(void *);
double ComputeForces( Particle [], Particle [], ParticleV [], int );
double ComputeNewPos( Particle [], ParticleV [], int, double);

Particle  * particles;   /* Particles */
ParticleV * pv;          /* Particle velocity */
int npart;

int main(int argc, char **argv)
{
    int         i, j;
    int         cnt;         /* number of times in loop */
    double      sim_t;       /* Simulation time */
    int tmp;
    if(argc != 3){
			printf("Wrong number of parameters.\nUsage: nbody num_bodies timesteps\n");
			exit(1);
		}

		npart = atoi(argv[1]);
		cnt = atoi(argv[2]);
		dt = 0.001;
		dt_old = 0.001;

    /* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);

    /* Generate the initial values */
    InitParticles(npart);
    sim_t = 0.0;

    while (cnt--) {
      double max_f;
      /* Compute forces (2D only) */
      max_f = ComputeForces( particles, particles, pv, npart );
      /* Once we have the forces, we compute the changes in position */
      sim_t += ComputeNewPos( particles, pv, npart, max_f);
    }
    for (i=0; i<npart; i++)
      fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);
    return 0;
}

PnV pnv;

void iniciar_particula(){
	pnv.ps.x	  = Random();
	pnv.ps.y	  = Random();
	pnv.ps.z	  = Random();
	pnv.ps.mass   = 1.0;
	pnv.psv.xold  = pnv.ps.x;
	pnv.psv.yold  = pnv.ps.y;
	pnv.psv.zold  = pnv.ps.z;
	pnv.psv.fx	  = 0;
	pnv.psv.fy	  = 0;
	pnv.psv.fz	  = 0;
}

void *thread_inicializadora(void *cont){
	int i = *((int *) cont);
	iniciar_particula();
	particles[i] = pnv.ps;
	pv[i] = pnv.psv;
	pthread_exit(NULL);
}

void InitParticles(int npart )
{
	pthread_t threads[npart];
    int i;
	int atrib[npart];
    for (i=0; i<npart; i++) {
		atrib[i] = i;
		pnv.ps = particles[i];
		pnv.psv = pv[i];
		pthread_create(&threads[i], NULL, thread_inicializadora, (void *)&atrib[i]);
    }
	for (i=0; i<npart; i++) {
		pthread_join(threads[i], NULL);
	}
}

double max_f;

void *thread_calculadora(void *cont){
    int i = *((int *) cont);
    int j;
    double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
    rmin = 100.0;
    xi   = particles[i].x;
    yi   = particles[i].y;
    fx   = 0.0;
    fy   = 0.0;
    for (j=0; j<npart; j++) {
      rx = xi - particles[j].x;
      ry = yi - particles[j].y;
      mj = particles[j].mass;
      r  = rx * rx + ry * ry;
      /* ignore overlap and same particle */
      if (r == 0.0) continue;
      if (r < rmin) rmin = r;
      r  = r * sqrt(r);
      fx -= mj * rx / r;
      fy -= mj * ry / r;
    }
    pv[i].fx += fx;
    pv[i].fy += fy;
    fx = sqrt(fx*fx + fy*fy)/rmin;
    if (fx > max_f) max_f = fx;
    pthread_exit(NULL);
}

double ComputeForces( Particle myparticles[], Particle others[], ParticleV pv[], int npart )
{
  int i;
  int atrib[npart];
  pthread_t threads[npart];
  int params[2];
  max_f = 0.0;
  for (i=0; i<npart; i++) {
      atrib[i] = i;
      pthread_create(&threads[i], NULL, thread_calculadora, (void *)&atrib[i]);
  }
  for (i=0; i<npart; i++) {
      pthread_join(threads[i], NULL);
  }
  return max_f;
}

void *thread_pos_calculadora(void *cont){
    int i = *((int *)cont);
    double a0, a1, a2;
    a0	 = 2.0 / (dt * (dt + dt_old));
    a2	 = 2.0 / (dt_old * (dt + dt_old));
    a1	 = -(a0 + a2);
    double xi, yi;
    xi	           = particles[i].x;
    yi	           = particles[i].y;
    particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
    particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
    pv[i].xold     = xi;
    pv[i].yold     = yi;
    pv[i].fx       = 0;
    pv[i].fy       = 0;
    pthread_exit(NULL);
}

double ComputeNewPos( Particle particles[], ParticleV pv[], int npart, double max_f)
{
  int i;
  int atrib[npart];
  double dt_new;
  pthread_t threads[npart];
  pthread_mutex_init(&mutex, NULL);
  for (i=0; i<npart; i++) {
      atrib[i] = i;
      pthread_create(&threads[i], NULL, thread_pos_calculadora, (void *)&atrib[i]);
  }
  for (i=0; i<npart; i++) {
      pthread_join(threads[i], NULL);
  }
  dt_new = 1.0/sqrt(max_f);
  /* Set a minimum: */
  if (dt_new < 1.0e-6) dt_new = 1.0e-6;
  /* Modify time step */
  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
    dt_old = dt;
    dt    *= 2.0;
  }
  pthread_mutex_destroy(&mutex);
  return dt_old;
}
