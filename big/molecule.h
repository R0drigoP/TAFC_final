#ifndef __molecule__
#define __molecule__

#include "TRandom3.h"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>

#include "nr3.h"
#include "mins.h"
#include "mins_ndim.h"

#include "global.h"
#include "config.h"

using namespace std;

class molecule{

public:

  //Constructors and Destructor
  molecule(unsigned int n_atoms = 0, float l_box = 0., float mute_prob = 0.);
  molecule(molecule* mom, molecule* dad, double gene_prop, unsigned int n_atoms = 0);
  ~molecule();

  //fitness
  void Fit();



  //Getters
  int Get_Natoms() {return N_atoms;}
  int Get_Dim() {return L_box;}
  double** Get_Pos() {return positions;}
  double Get_Fit() {return fitness;}

  //Setters
  void Set_Pos(double** pos);
  void Set_Fit(double pot);
  
  //Comparator
  static bool LessPot(molecule* mol1, molecule* mol2) {return mol1->Get_Fit() < mol2->Get_Fit();}

  //Mutation and Reproduction
  void Mutate(unsigned int iter, float m0, float alpha, int flag, TRandom3* gRandom = NULL);
  int generate_children(vector<molecule*> pop, TRandom3* gRandom);
  void Mating(molecule* mom, molecule* dad, double gene_prop = 0.);
  void Mating_Plano(molecule* mom, molecule* dad, TRandom3* gRandom);
  void Mating_Sphere(molecule* mom, molecule* dad, TRandom3* gRandom);



private:
  unsigned int N_atoms;
  float mutation_prob, L_box;
  double **positions;
  double fitness;

  void Local_Min();
};

#endif

struct Funcd {
  
    Funcd(){}

    Doub operator() (VecDoub_I &x) {
      Doub f = 0.;
      for (int i = 0; i < N_atoms-1; ++i){
          for(int j = i+1; j < N_atoms; ++j){
              // Calculate radius
              Doub r = 0.;
              for(int k = 0; k < 3; ++k)
                  r += (x[i*3+k] - x[j*3+k])*(x[i*3+k] - x[j*3+k]);
              Doub inv =1/r;
              Doub inv3 = inv*inv*inv;

              f += 4*(inv3*inv3- inv3);
          }
      }
      nb_of_calls ++;
      return f;
    }

    void df(VecDoub_I &x, VecDoub_O &deriv) {

        for (int i = 0; i < N_atoms-1; ++i){
            for(int j = i+1; j < N_atoms; ++j){
                // Calculate radius
                Doub r = 0.;
                for(int k = 0; k < 3; ++k)
                    r += (x[i*3+k] - x[j*3+k])*(x[i*3+k] - x[j*3+k]);
                Doub inv = 1/r;
                Doub sqr = 1/sqrt(r);
                Doub inv3 = inv*inv*inv;


                // Calculate force and add to deriv
                Doub force = -24*sqr*inv3*( 2*inv3   - 1 );
                for(int k = 0; k < 3; ++k){
                    deriv[i*3+k] += force*(x[i*3+k] - x[j*3+k]);
                    deriv[j*3+k] += force*(x[j*3+k] - x[i*3+k]);
                }
            }
        }
        nb_of_calls_der ++; 
    }
};
