#include "molecule.h"


//Constructors and Destructor
molecule::molecule(unsigned int n_atoms, float l_box, float mute_prob) : N_atoms(n_atoms), L_box(l_box), mutation_prob(mute_prob) {
  positions = new double*[n_atoms];
  
  for (int i = 0; i < n_atoms; i++)
    positions[i] = new double[3];
  
  
  //random positions in a box
  if(initial_structure == 0){
    
    gRandom = new TRandom3(0);
    
    for (int i = 0; i < N_atoms; i++){
      for (int j = 0; j < 3; j++)
	positions[i][j] = gRandom -> Uniform(0., l_box);
    }
    delete gRandom;
  }
  
  //spherical surface inside another spherical surface + atom in centre
  if(initial_structure == 1){
    
    //should be ajusted for different cluster sizes
    int N_atoms_out = int(N_atoms*0.7);

    double r_out = l_box/2.;
    double r_in =  l_box/4.;

    //atom in centre
    positions[0][0] = 0.;
    positions[0][1] = 0.;
    positions[0][2] = 0.;

    for (int i = 1; i < N_atoms_out; ++i){

      double theta = gRandom -> Uniform(0., 2*M_PI);
      double phi = gRandom -> Uniform(0., M_PI);
      

      positions[i][0] = r_out*sin(theta)*cos(phi);
      positions[i][1] = r_out*sin(theta)*sin(phi);
      positions[i][2] = r_out*cos(theta);
    }

    for (int i = N_atoms_out; i < N_atoms; ++i){
      double theta = gRandom -> Uniform(0., 2*M_PI);
      double phi = gRandom -> Uniform(0., M_PI);

      positions[i][0] = r_in*sin(theta)*cos(phi);
      positions[i][1] = r_in*sin(theta)*sin(phi);
      positions[i][2] = r_in*cos(theta);
    }
  }

  this->Fit();

}

molecule::molecule(molecule* mom, molecule* dad, double gene_prop, unsigned int n_atoms): N_atoms(n_atoms) {
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();
  
  for (int i = 0; i < N_atoms; i++)
    for (int j = 0; j < 3; j++)
      positions[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
    
}

molecule::~molecule() {

  for(int i = 0; i < 3; ++i){ 
    delete[] positions[i];
  }

  delete[] positions;
}


//Method that calculates fitness
void molecule::Fit(){
  double f = 0.;

#pragma omp atomic
  ++nb_of_calls;

  for(int i = 0; i < N_atoms-1; ++i){
    for(int j = i+1; j < N_atoms; ++j){

      //Calculate radius
      double r_square = 0.;

      for(int k = 0; k < 3; ++k)
        r_square += (positions[i][k] - positions[j][k])*(positions[i][k] - positions[j][k]);

      double r_four = r_square*r_square;

      //Calculate Potential and sum to f_value
      f += 4*( 1/(r_four*r_four*r_four) - 1/(r_four*r_square));
    }
  }

  fitness = f;
}

//Setters
void molecule::Set_Pos(double** new_pos){
  for(int i = 0; i < N_atoms; ++i){
    for(int j = 0; j < 3; ++j)
      positions[i][j] = new_pos[i][j];
  }
}

void molecule::Set_Fit(double fit){
  fitness = fit;
}

//Mutation and Reproduction
void molecule::Mutate(unsigned int iter, float m0, float alpha, int flag, TRandom3* gRandom){
  
  double check_if_mute = gRandom->Uniform(0,1);

  if(check_if_mute < mutation_prob){

    int atom_to_mutate = (int)gRandom->Uniform(0, N_atoms);

    m0 = m0/(1 + alpha*iter);

    if(m0 < 0.000001)
      m0 = 0.000001;

    for(int i = 0; i < 3; ++i){

      double x0 = gRandom->Uniform(-1,1)*m0*L_box;

      positions[atom_to_mutate][i] += x0;
    
    }
    this->Fit();
  }

  //In case the individual didnt mutate it calcules the potential for the ones that had sexual reproduction
  else if(flag == 1)
    this->Fit();
}

int molecule::generate_children(vector<molecule*> pop, TRandom3* gRandom){

  double check_if_sex = gRandom->Uniform(0,1);

  if(check_if_sex < sex_prob ){

    //choose random parents from the surviving population
    int mom_index = (int)gRandom->Uniform(0, parents_nb);
    int dad_index = (int)gRandom->Uniform(0, parents_nb);

    //and make sure they are different
    while(dad_index==mom_index)
      dad_index = (int)gRandom->Uniform(0, parents_nb);

    double mating_type = gRandom->Uniform(0,1);

    //call reproduction method
    if(mating_type < 0.9){
      this -> Mating_Plane(pop[mom_index],pop[dad_index], gRandom);
    }

    else{
      this -> Mating(pop[mom_index],pop[dad_index],gRandom->Uniform(0,1));
    }
    
    return 1;
  }
  
  return 0;  
}

void molecule::Mating(molecule* mom, molecule* dad, double gene_prop){
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();
  
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++){
      positions[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
    }
  }
}

void molecule::Mating_Plane(molecule* mom, molecule* dad, TRandom3* gRandom){

  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();
  double pos_CM = 0;

  int flag = 0;
  int nr_atoms = 0;

  int dir = (int)gRandom -> Uniform(0,3);

  //centre of mass 
  for(int i = 0; i < N_atoms; i++)
    pos_CM += (pos_mom[i][dir] + pos_dad[i][dir])/2;
  pos_CM = pos_CM/N_atoms; 

  //choose every mom atom above CM
  for (int i = 0; i < N_atoms; i++){

    if(pos_mom[i][dir] > pos_CM){
      for (int k = 0; k < nr_atoms; k++){
        flag = 0;
        for (int j = 0; j < 3; j++){
          if (positions[k][j]==pos_mom[i][j])
            flag++;
        }
        if(flag == 3)
          break; 
      }
      if (flag != 3){
	for(int j = 0; j < 3; j++)
	  positions[nr_atoms][j] = pos_mom[i][j];
	nr_atoms ++;
	//cout<<"mom  above"<<endl;
      }
    } 
  }

  //choose every dad atom bellow CM
  for (int i = 0; i < N_atoms; i++){
    if(nr_atoms== N_atoms)
      break;
    else if(pos_dad[i][dir] < pos_CM){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_dad[i][j];
      nr_atoms ++;
    }
  }

  if (nr_atoms < N_atoms){
    //choose every mom atom below CM
    for (int i = 0; i < N_atoms; i++){

      if(pos_mom[i][dir] < pos_CM){
        for (int k = 0; k < nr_atoms; k++){
          flag = 0;
          for (int j = 0; j < 3; j++){
            if (positions[k][j]==pos_mom[i][j])
              flag++;
          }
          if(flag == 3)
            break; 
        }
        if (flag != 3){
          for(int j = 0; j < 3; j++)
            positions[nr_atoms][j] = pos_mom[i][j];
          nr_atoms ++;
          //cout<<"mom below"<<endl;
        }
      }
      if(nr_atoms== N_atoms)
        break;
    }
  }

  if (nr_atoms < N_atoms){
    //choose every dad atom above CM
    for (int i = 0; i < N_atoms; i++){

      if(pos_dad[i][dir] > pos_CM){
        for (int k = 0; k < nr_atoms; k++){
          flag = 0;
          for (int j = 0; j < 3; j++){
            if (positions[k][j]==pos_dad[i][j])
              flag++;
          }
          if(flag == 3)
            break; 
        }
        if (flag != 3){
          for(int j = 0; j < 3; j++)
            positions[nr_atoms][j] = pos_dad[i][j];
          nr_atoms ++;
          //cout<<"dad above"<<endl;
        }
      }
      if(nr_atoms == N_atoms)
        break;
    }
  }
}
