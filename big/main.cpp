

#include "molecule.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include <ctime>
#include <fstream>
#include "parameters.h"
//#include <omp.h>


//

using namespace std;

//global variables

unsigned int parents_nb = int(survival_rate * N_molecules);

unsigned int nb_of_calls = 0, nb_of_calls_der=0, nb_of_calls_mute = 0, nb_of_calls_mat = 0, nb_of_calls_mat_plano = 0;

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

  if(N_atoms_out+N_atoms_in != N_atoms){
    cout<<"Please specify a coherent number of atoms in config.h"<<endl;
    return 1;
  }

  if (local_min_iter > max_iter){
    cout<<"Please be more coherent about at which iteration you start the local minimization after mutations"<<endl;
    return 1;
  }

  double **positions;
  positions = new double*[N_atoms];
  for (int i = 0; i < N_atoms; i++) 
    positions[i] = new double[3];

  int* flag = new int[N_molecules];
  for (int i = 0; i < N_molecules; i++) 
    flag[i] = 0;


  //TCanvas *c1 = new TCanvas();
  //auto gr = new TGraph();

  ofstream text_file("best_molecule.bs");
  ofstream movie_file("best_molecule.mv");

   ofstream opt_file("opt.txt");//used for communication with optuna

  
  //population of molecules
  vector<molecule*> pop(N_molecules);

  for(int i = 0; i < N_molecules; ++i)
    pop[i] = new molecule(N_atoms, L_box, mutation_prob);

  Funcd funcd;
  Frprmn<Funcd> frprmn(funcd);
  VecDoub p(3*N_atoms);

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
          flag[i] = pop[i]->generate_children(pop, gRandom);

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

    

    //matar os mais fracos e fazer copias da melhor pop (se calhar atribuir alguma aleatoriadade a este processo)

      for(int mol = static_cast<int>(survival_rate*N_molecules); mol < N_molecules; mol += static_cast<int>(survival_rate*N_molecules)){
        int alive = 0;
        while(alive < survival_rate*N_molecules && (mol+alive)<N_molecules){
        //cout<<mol<<" "<<alive<<endl;                          
          pop[mol+alive] -> Set_Pos(pop[alive] -> Get_Pos());
          pop[mol+alive] -> Set_Fit(pop[alive] -> Get_Fit());

          ++alive;
        }
      }
    

    //print to movie file
    
    
    if(iter == 0){
      double** best_pos = pop[0]->Get_Pos();
      for(int i = 0; i < N_atoms; ++i){
        text_file << "atom C " << flush;
        for(int j = 0; j < 3; ++j)
          text_file << best_pos[i][j] << " " << flush;
        text_file << endl;
      }
    }
    else if(iter % 100 == 0){
      movie_file << "frame" << endl;
      double** best_pos = pop[0]->Get_Pos();
      for(int i = 0; i < N_atoms; ++i){
        for(int j = 0; j < 3; ++j)
          movie_file << best_pos[i][j] << " ";
      }
      movie_file << endl << endl;
    }
    
  }//closing iterations loop


  cout << "Pot: " << pop[0] -> Get_Fit() << endl;

  opt_file<<pop[0] -> Get_Fit()<<flush;
  opt_file.close();

  
  for(int i = 0; i < 3; ++i) 
    delete[] positions[i];
  
  delete[] positions;
  delete[] flag;
  
  pop.clear();

  
  double atom_size = 0.1/L_box;
  
  text_file << endl << "spec C 0.1 Red" << endl
	    << endl << "bonds C C 0.3 1.5 0.01 0.0"
	    << endl << "bonds C H 0.4 1.0 0.01 1.0"
	    << endl << "scale 100"
	    << endl << "inc 5.0" <<endl<< flush;
  
/*
  c1 -> cd();
  gr->GetHistogram()->SetMaximum(-0.1*final_fit);
  gr->GetHistogram()->SetMinimum(1.01*final_fit);

  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");*/
  
  text_file.close();
  movie_file.close();


  cout<<"Total pot calls: "<<nb_of_calls<<endl;
  cout<<"Total der calls: "<<nb_of_calls_der<<endl;


  /*
  double t1 = omp_get_wtime();
  printf("\n");
  printf("Number of threads   =  %i\n", omp_get_max_threads());
  printf("Computation time    =  %f ms\n", (t1-t0) * 1000);  */
  return 0;
}
