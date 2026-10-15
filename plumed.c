#include "elphy.h"

#ifdef _PLUMED
static int nph;
static double *pos, *pos0;
static plumed p;

void plumed_init(int *nat, double *dt, double (*uc)[3], double (*tau)[3],
    char *dat, char *log) {

    static int prec = sizeof(double);
    static double one = 1.0;

    p = plumed_create();

    plumed_cmd(p, "setRealPrecision", &prec);

    plumed_cmd(p, "setMDEnergyUnits", &one);
    plumed_cmd(p, "setMDLengthUnits", &one);
    plumed_cmd(p, "setMDTimeUnits", &one);

    plumed_cmd(p, "setNatoms", nat);
    plumed_cmd(p, "setTimestep", dt);
    plumed_cmd(p, "setNoVirial", NULL);

    plumed_cmd(p, "setMDEngine", "elphy");
    plumed_cmd(p, "setPlumedDat", dat);
    if (log)
        plumed_cmd(p, "setLogFile", log);

    plumed_cmd(p, "init", NULL);

    plumed_cmd(p, "setBox", *uc);

    nph = 3 * *nat;
    pos0 = *tau;

    if (!(pos = malloc(nph * sizeof *pos)))
        error("No memory for atomic positions.");
}

void plumed_bias(const double *u, double *forces, double *energy) {
    const int inc = 1;
    const double plus = 1.0;
    static int i = 0;

    plumed_cmd(p, "setStep", &i);
    i++;

    memcpy(pos, pos0, nph * sizeof *pos);
    daxpy_(&nph, &plus, u, &inc, pos, &inc);

    plumed_cmd(p, "setPositions", pos);
    plumed_cmd(p, "setEnergy", energy);
    plumed_cmd(p, "setForces", forces);

    plumed_cmd(p, "calc", NULL);
}

void plumed_exit(void) {
    plumed_finalize(p);
    free(pos);
}
#endif
