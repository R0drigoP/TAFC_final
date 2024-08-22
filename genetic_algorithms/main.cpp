#include "molecule.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include <ctime>
#include <fstream>
#include "global.h"

using namespace std;

unsigned int N_molecules = 30, N_atoms = 13;
float L_box = .8, survival_rate = 0.9;
unsigned int max_iter = 100000;

//mutation parameters
float mutation_prob = 0.01, alpha = 0, m0 = 0.01;

//mating, if sex_prob = 0, no mating will occur
bool mating = 0;
float sex_prob = 0.;

//inicial_structure = 0 : Random positions in box
//inicial_structure = 1 : Random positions in 2 spherial surface (one inside the other and 1 atom in centre)
bool initial_structure = 0;

unsigned int parents_nb = int(survival_rate * N_molecules);

//nr of iterations with similar best fitness before we stop our algorithm
int acceptance = 35000;

//In order to obtain the potential with six decimal places local_minization should be 1 
//However with this minimization it becomes more difficult to arrive at the global mininum
//since it is likely to be stuck in local minimums.
bool minimization = 0;
unsigned int iter_same_potential = 0;
unsigned int nb_of_calls = 0;

int main(){

  double final_fit = 0.;
  int minimum = 0;

  //best fitness obtained 
  float best;

  double **positions;
  positions = new double*[N_atoms];
  for (int i = 0; i < N_atoms; i++) 
    positions[i] = new double[3];

  int* flag = new int[N_molecules];
  for (int i = 0; i < N_molecules; i++) 
    flag[i] = 0;

  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();

  ofstream text_file("best_molecule.bs");
  ofstream movie_file("best_molecule.mv");
  
  //population of molecules
  vector<molecule*> pop(N_molecules);

  for(int i = 0; i < N_molecules; ++i)
    pop[i] = new molecule(N_atoms, L_box, mutation_prob);

  best = pop[0] -> Get_Fit();

  for(int iter = 0; iter < max_iter; iter++){

    if(best > pop[0] -> Get_Fit())
      best = pop[0] -> Get_Fit();

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


    //kill weaker individuals e make copies of the best individuals
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
    else if(iter % 1000 == 0){
      movie_file << "frame" << endl;
      double** best_pos = pop[0]->Get_Pos();
      for(int i = 0; i < N_atoms; ++i){
        for(int j = 0; j < 3; ++j)
          movie_file << best_pos[i][j] << " ";
      }
      movie_file << endl << endl;
    }

    gr -> AddPoint( iter, pop[0] -> Get_Fit());

    //prints
    /*if(iter%20000 == 0){
      cout << "Iteration nr : " << iter << endl;
      cout << "Potencial :  " << setprecision(8) << pop[0] -> Get_Fit() << endl;
    }

    if(iter == max_iter-1)
      cout << "Final Potential : " << setprecision(8) << pop[0]-> Get_Fit() << endl;
    
    */

    //early stopping
    if( (pop[0] -> Get_Fit()) < best)
      iter_same_potential = 0;
    else{
      iter_same_potential++;
      if( iter_same_potential == acceptance)
        break;
    }

    if(minimization == 1){
      if(((best - pop[0] -> Get_Fit()) < m0/1000)){
        minimum++;
      }

      if(((best - pop[0] -> Get_Fit()) > m0/1000)){
        minimum = 0;
      }

      if(minimum > 10000 && m0 > 0.0000001){
        mutation_prob = mutation_prob*2;
        mating = 0;
        m0 = m0/10;
        minimum = 0;
      }
    }

  }//closing iterations loop

  cout << "Final fitness: " << pop[0] -> Get_Fit() << endl;

  for(int i = 0; i < 3; ++i) 
    delete[] positions[i];
  
  delete[] positions;
  delete[] flag;

  pop.clear();
  
  text_file << endl << "spec C 0.1 Red" << endl
	    << endl << "bonds C C 0.3 1.5 0.01 0.0"
	    << endl << "bonds C H 0.4 1.0 0.01 1.0"
	    << endl << "scale 100"
	    << endl << "inc 5.0" << endl;
  

  c1 -> cd();
  c1 -> SetLogx();

  gr -> SetMarkerStyle(8);
  gr -> SetMarkerColor(kRed + 1);
  gr -> SetMarkerSize(0.3);
    
  gr -> GetYaxis() -> SetTitle("Potential");
  gr -> GetXaxis() -> SetTitle("Function Evaluations");

  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");
  
  text_file.close();
  movie_file.close();

  delete gr;
  delete c1;

  cout << "Potential calls: " << nb_of_calls << endl;

  return 0;
}
