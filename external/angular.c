/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <time.h>


double computeCutOff(double rij, double rcut) {
    /* computeCutOff */
    if (rij>rcut){
        return 0.; }
    if (rij<=rcut){
        return 0.5*(cos((rij*3.14159265359)/rcut)+1); }
}

double computeG4(double sqrij, double sqrik, double sqrjk, double eta, double fcij, double fcik, double fcjk ) {
    /* computeG4 */
    double fc4;
    fc4 = fcij * fcik * fcjk;
    return exp( -eta*( sqrij+sqrik+sqrjk ) ) * fc4; }

double computeG5(double sqrij, double sqrik, double eta, double fcij, double fcik ) {
    /* computeG5 */
    double fc5;
    fc5 = fcij * fcik;
    return  exp( -eta*( sqrij+sqrik ) ) * fc5; }

int sindex(int *spe, int val) {
    int j;
    for(j=0;j<sizeof(spe);j++){ 
        if (spe[j]==val){return j;}
    }
}

void computeAngular(
        int natom,
        int index_i,
        int nspecies,
        int nAngular,
        double  rcut,
        double *eta_g4, 
        double *zet_g4, 
        double *lbd_g4,
        double *eta_g5, 
        double *zet_g5, 
        double *lbd_g5,
        int    *species,
        int    *numSymbols, 
        double *distanceMatrix,
        double *g4,
        double *g5 ) {
    /* computeAngular */

    int j,k,h;
    double hZj,tmp;

    double rij,rik,rjk;
    double fcij,fcik,fcjk;
    int index_j,index_k;

    double costheta;
    double gauss_g4,gauss_g5;

    int n = nspecies ;
    int l = nAngular ;

    for(j=0;j<natom;j++){ 
        for(k=0;k<natom;k++){
            if(index_i==k){ continue; }
            if(index_i==j){ continue; }
            if(k>=j){ continue; }

            rij = distanceMatrix[index_i*natom + j];
            
            if(rij==0.0){continue;}
            if(rij>rcut){continue;}


            //hZj     = numSymbols[j] * numSymbols[k];
            hZj     = numSymbols[j] * numSymbols[k] / (numSymbols[j] + numSymbols[k]);
            
            index_j = sindex(species, numSymbols[j]) ;
            index_k = sindex(species, numSymbols[k]) ;

            //printf("%d %d \n",index_j,index_k);

            rik = distanceMatrix[index_i*natom + k];
            rjk = distanceMatrix[j*natom + k];

            fcij = computeCutOff(rij ,rcut );
            fcik = computeCutOff(rik ,rcut );
            fcjk = computeCutOff(rjk ,rcut );

            costheta = 0.5/(rij*rik) * (pow(rij,2) + pow(rik,2) - pow(rjk,2) );

            for(h=0;h<l;h++) {
                gauss_g4  = computeG4( pow(rij,2), pow(rik,2), pow(rjk,2), eta_g4[h], fcij, fcik, fcjk) ;
                tmp       = 0.5*(1 + lbd_g4[h] * costheta) ;
                g4[index_j*(n*l)+index_k*(l)+h] += 2 * pow(tmp, zet_g4[h]) * gauss_g4 * hZj;
                if(index_j!=index_k) {
                    g4[index_k*(n*l)+index_j*(l)+h] += 2 * pow(tmp, zet_g4[h]) * gauss_g4 * hZj;
                }
            
                gauss_g5  = computeG5( pow(rij,2), pow(rik,2), eta_g5[h], fcij, fcik) ;
                tmp       = 0.5*(1 + lbd_g5[h] * costheta) ;
                g5[index_j*(n*l)+index_k*(l)+h] += 2 * pow(tmp, zet_g5[h]) * gauss_g5 * hZj;
                if(index_j!=index_k) {
                    g5[index_k*(n*l)+index_j*(l)+h] += 2 * pow(tmp, zet_g5[h]) * gauss_g5 * hZj;
                }
            }
        }
    }
}

