
#include "molecule.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include <ctime>
#include <fstream>



using namespace std;

//global variables
unsigned int N_atoms = 13;




unsigned int parents_nb = int(survival_rate * N_molecules);

unsigned int nb_of_calls = 0, nb_of_calls_mute = 0, nb_of_calls_mat = 0, nb_of_calls_mat_plano = 0;

double final_fit = 0.;
int main(){

  if(mating==1 &&  N_molecules * survival_rate < 2.){
    cout<<"To have sexual reprodution at least 2 molecules must survive each gen..."<<endl;
    cout<<"Increase your population or the survival probability"<<endl;
    return 1;
  }

  if(survival_rate < 0. || survival_rate > 1. || parents_nb >= N_molecules){
    cout << "Check your survival_rate" << endl;
    return 1;
  }

  double **positions;
  positions = new double*[N_atoms];
  for (int i = 0; i < N_atoms; i++) 
    positions[i] = new double[3];

  int* flag = new int[N_molecules];
  for (int i = 0; i < N_molecules; i++) 
    flag[i] = 0;

  ofstream opt_file("opt.txt");
  ofstream opt_calls("calls.txt");

  
  //population of molecules
  vector<molecule*> pop(N_molecules);

  for(int i = 0; i < N_molecules; ++i)
    pop[i] = new molecule(N_atoms, L_box, mutation_prob);

  Funcd funcd;
  Frprmn<Funcd> frprmn(funcd);
  VecDoub p(3*N_atoms);

  //variables for early stopping
  double best = 0.; //so far best potential
  int iter_stop = 0; //iterations with the same best so far
  int acceptance = 1000; //nb of iterations allowed in the same best

  for(int iter = 0; iter < max_iter; iter++){


    //REPRODUCTION 

    #pragma omp parallel
    {
      //---sexual reproduction
      if( mating == 1){

       
        for(int i = 0 ; i < parents_nb; i++)
          flag[i] = 0; //setting flag to 0 for parents

        TRandom3* gRandom = new TRandom3(0); 


        #pragma omp for
        for(int i = parents_nb ; i <N_molecules; i++)
          flag[i] = pop[i]->generate_children3(pop, gRandom);

        delete gRandom;
      }

      //---assexual reproduction

      TRandom3* gRandom = new TRandom3(0); 
      #pragma omp for
      for(int mol = 0; mol < N_molecules; mol++)
        pop[mol] -> Mutate(iter, m0, alpha, flag[mol], gRandom);
      
      delete gRandom;

    }

    //sort population
    sort(pop.begin(), pop.end(), molecule::LessPot);

    //matar os mais fracos e fazer copias da melhor pop
    #pragma omp parallel
    {
      #pragma omp for
      for(int mol = static_cast<int>(survival_rate*N_molecules); mol < N_molecules; mol += static_cast<int>(survival_rate*N_molecules)){
        int alive = 0;
        while(alive < survival_rate*N_molecules && (mol+alive)<N_molecules){            
          pop[mol+alive] -> Set_Pos(pop[alive] -> Get_Pos());
          pop[mol+alive] -> Set_Fit(pop[alive] -> Get_Fit());

          ++alive;
        }
      }
    }

    //early stopping
    double current = pop[0]->Get_Fit();
    if (iter==0)
      best = current;
    if ( current < best){
      iter_stop=0;
      best = current;
    }
    else{
      iter_stop ++;
      if(iter_stop == acceptance  )
        break;
    }
    
  }//closing iterations loop


  cout << "Final fitness: " << pop[0] -> Get_Fit() << endl;

  opt_file<<pop[0] -> Get_Fit()<<flush;
  opt_file.close();

  opt_calls<<nb_of_calls<<flush;
  opt_calls.close();
  
  for(int i = 0; i < 3; ++i) 
    delete[] positions[i];
  
  delete[] positions;
  delete[] flag;
  
  pop.clear();
  

  cout<<"Total pot calls: "<<nb_of_calls<<endl;
  return 0;
}
