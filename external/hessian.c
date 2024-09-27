#include <stdio.h>
#include <math.h>
#include <time.h>
//#include "cblas"

void mkde2(int natom, int ngeff, int *ngindx, double *wn, double *mff, double *mwq, double *gg, double *de2)
{
   int mg,ng;
   int a,b,x,y;
   int iat,jat,kat,hat;
   int k,l,h,g,f,d,e;
   double u,v,w;


   clock_t begin;
   clock_t end;
   double sec = 0.;
   double seconde = 0.;

   double wc[natom*natom];

   double wa[natom];
   double wb[natom];
   
   int n_idx  = ngeff * 2 ;
   int n_atom = natom * natom;

   /*
   double mw[natom*natom];
   
   for(iat=0;iat<ngeff*natom;iat++) {
      for(jat=0;jat<ngeff*natom;jat++) {
         mw[iat%natom+jat%natom] += mwq[iat] * mwqt[jat];
      }
   }
   
   //extern void dgemm();
   //dgemm(TRANSA="N", TRANSB="N", A=mwq, B=mwq, C=pHess );  

   for(iat=0;iat<ngeff*natom;iat++) {
      printf("%e ",mw[iat]);
   }*/
   

   for(mg=0;mg<n_idx;mg+=2){
      
      for(ng=mg;ng<n_idx;ng+=2){

         a = ngindx[mg];
         b = ngindx[ng];
         x = ngindx[mg+1];
         y = ngindx[ng+1];

         u = 0.0;
         if(a == b) u += mff[x*natom+y]*( wn[a]+wn[b]-wn[x]-wn[y]);
         if(x == y) u += mff[a*natom+b]*(-wn[a]-wn[b]+wn[x]+wn[y]);
         if(x == b) u -= mff[a*natom+y]*(-wn[a]+wn[b]+wn[x]-wn[y]);
         if(a == y) u -= mff[b*natom+x]*( wn[a]-wn[b]-wn[x]+wn[y]);

         v = 0.0;
   
         /*
         for(iat=0;iat<natom;iat++) {
            g = iat*n_atom+e;
            h = iat*n_atom+d;
            wa[iat] = mwq[g];
            wb[iat] = mwq[h];
         }

         for(iat=0;iat<natom;iat++) {
            for(jat=0;jat<natom;jat++) {
               f = iat*natom+jat;
               wc[f] = wa[iat] * wb[jat];
            }
         }

         for(iat=0;iat<natom*natom;iat++) {
            v += wc[iat] * gg[iat] ;   
         }
         */
         
         d = a*natom+x;
         e = b*natom+y;
         for(iat=0;iat<natom;iat++) {
               g = iat*n_atom+e;
               w = mwq[g];
               for(jat=0;jat<natom;jat++) {
                  h = jat*n_atom+d;
                  f = iat*natom+jat;
                  v += w*gg[f]*mwq[h];
               }
            }
         
         k = mg/2;
         l = ng/2;

         v *= 4.0*(wn[a] - wn[x])*(wn[b] - wn[y]);

         de2[k*ngeff + l] = u + v;
         de2[l*ngeff + k] = u + v;
         

      }
   }
}


