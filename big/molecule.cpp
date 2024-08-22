#include "molecule.h"


//Constructors and Destructor

molecule::molecule(unsigned int n_atoms, float l_box, float mute_prob) : N_atoms(n_atoms), L_box(l_box), mutation_prob(mute_prob) {
  positions = new double*[n_atoms];
  
  for (int i = 0; i < n_atoms; i++)
    positions[i] = new double[3];
  
  gRandom = new TRandom3(0);

  /* random initialization in 1 spherical surface
  for (int i = 0; i < N_atoms; ++i){
      //double r = gRandom -> Uniform(0., l_box/2.);
      double r = l_box/2.;
      double theta = gRandom -> Uniform(0., 2*M_PI);
      double phi = gRandom -> Uniform(0., M_PI);
      
      positions[i][0] = r*sin(theta)*cos(phi);
      positions[i][1] = r*sin(theta)*sin(phi);
      positions[i][2] = r*cos(theta);
    }*/

  /* random initialization in a box
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++)
    positions[i][j] = gRandom -> Uniform(0., l_box);
  }*/

  
  //random initialization in 2 spherical concentric spheres


      double r_in  = l_box/3.;
      double r_out =  3*l_box/4.;

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

  this->Fit();
}

molecule::molecule(molecule* mom, molecule* dad, double gene_prop, unsigned int n_atoms): N_atoms(n_atoms) {
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();
  
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++){
      positions[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
    }
  }
}

molecule::~molecule() {
  for(int i = 0; i < 3; ++i){ 
    delete[] positions[i];
  }
  
  delete[] positions;
}

//Private method that calculates fitness
void molecule::Fit(){
  double f = 0.;

  ++nb_of_calls;

  for(int i = 0; i < N_atoms-1; ++i){
    for(int j = i+1; j < N_atoms; ++j){
      //Calculate radius
      double r = 0.;
      for(int k = 0; k < 3; ++k)
        r += (positions[i][k] - positions[j][k])*(positions[i][k] - positions[j][k]);

      double inv = 1./r;
      double inv3=inv*inv*inv;

      //Calculate Potential and sum to f_value
      f += 4*inv3*( inv3 - 1);
    }
  }
  fitness = f;
}

//Setters

void molecule::Set_Pos(double** new_pos){
  for(int i=0; i<N_atoms; ++i){
    for(int j=0; j<3; ++j)
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
    nb_of_calls_mute ++;
    int atom_to_mutate = (int)gRandom->Uniform(0, N_atoms);

    for(int i=0; i<3; ++i){
      
      double x0 = gRandom->Uniform(-1,1)*m0*L_box;

      double mutation = x0/(1+alpha*iter);

      positions[atom_to_mutate][i] += mutation;
      
    }
    if(iter> local_min_iter)
      this->Local_Min();


    this->Fit();
  }
  //se nao tiver feito mutacao calcula na mesma o potencial para os que sofreram reproducao sexuada
  else if(flag==1)
    this->Fit();
}

void molecule::Local_Min(){

  VecDoub ploc(3*N_atoms);
  Funcd funcd;
  Frprmn<Funcd> frprmn(funcd);

  for(int i = 0; i < N_atoms; i++ ){
    for(int j = 0; j < 3; j++){
      ploc[i*3+j] = positions[i][j];
    }
  }

  ploc = frprmn.minimize(ploc);

  for(int i = 0; i < N_atoms; i++ ){
    for(int j = 0; j < 3; j++){
      positions[i][j] = ploc[i*3+j];
    }
  }

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
    if(mating_type<0.9){
      nb_of_calls_mat_plano ++;
      this->Mating_Plano(pop[mom_index],pop[dad_index], gRandom);
    }

    else{
      nb_of_calls_mat ++;
      this->Mating(pop[mom_index],pop[dad_index],gRandom->Uniform(0,1));
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

void molecule::Mating_Sphere(molecule* mom, molecule* dad, TRandom3* gRandom){

  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();

  //min number of atoms of each parent
  const int minAtoms = 2;
  int nAtomsMom = 0;
  int nAtomsDad = 0;
  double r = 0., x=0., y=0., z=0.;

  int flag = 0;
  int nr_atoms = 0;

  //setting positions to 0
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++)
    positions[i][j] = 0;
  }

  //checking that there are at least 2 atoms for each parent inside the sphere
  while (nAtomsMom <= minAtoms || nAtomsDad <= minAtoms) {
    //randomly selecting the radius 
    r = gRandom->Uniform(0,L_box);

    //center coordinates of sphere using normal dist
    x = gRandom->Gaus();
    y = gRandom->Gaus();
    z = gRandom->Gaus();

    //reset to 0
    nAtomsMom = 0;
    nAtomsDad = 0;

    //loop over atoms
    for (int i = 0; i < N_atoms; i++){
      double dist_mom = sqrt((pos_mom[i][0] - x)*(pos_mom[i][0] - x) + (pos_mom[i][1] - y)*(pos_mom[i][1] - y) + (pos_mom[i][1] - z)*(pos_mom[i][1] - z));
      double dist_dad = sqrt((pos_dad[i][0] - x)*(pos_dad[i][0] - x) + (pos_dad[i][1] - y)*(pos_dad[i][1] - y) + (pos_dad[i][1] - z)*(pos_dad[i][1] - z));
      if (dist_mom <= r) 
        nAtomsMom++;
      if (dist_dad <= r) 
        nAtomsDad++;
    }
  }

  //choose every mom atom inside sphere
  for (int i = 0; i < N_atoms; i++){

    double dist_mom = sqrt((pos_mom[i][0] - x)*(pos_mom[i][0] - x) + (pos_mom[i][1] - y)*(pos_mom[i][1] - y) + (pos_mom[i][1] - z)*(pos_mom[i][1] - z));

    if(dist_mom <= r){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_mom[i][j];
      nr_atoms ++;
    }
  }

  //choose every dad atom outside sphere
  for (int i = 0; i < N_atoms; i++){
    double dist_dad = sqrt((pos_dad[i][0] - x)*(pos_dad[i][0] - x) + (pos_dad[i][1] - y)*(pos_dad[i][1] - y) + (pos_dad[i][1] - z)*(pos_dad[i][1] - z));
    if(nr_atoms== N_atoms)
      break;
    else if( dist_dad > r){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_dad[i][j];
      nr_atoms ++;
    }
  }
  if (nr_atoms < N_atoms){
    //choose every mom outside 
    for (int i = 0; i < N_atoms; i++){
      double dist_mom = sqrt((pos_mom[i][0] - x)*(pos_mom[i][0] - x) + (pos_mom[i][1] - y)*(pos_mom[i][1] - y) + (pos_mom[i][1] - z)*(pos_mom[i][1] - z));

      if( dist_mom > r){
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
        }
      }
      if(nr_atoms== N_atoms)
        break;
    }
  }

  if (nr_atoms < N_atoms){
    //choose every dad atom above CM
    for (int i = 0; i < N_atoms; i++){
      double dist_dad = sqrt((pos_dad[i][0] - x)*(pos_dad[i][0] - x) + (pos_dad[i][1] - y)*(pos_dad[i][1] - y) + (pos_dad[i][1] - z)*(pos_dad[i][1] - z));

      if(dist_dad <= r){
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
      if(nr_atoms== N_atoms)
        break;
    }
  }
}




void molecule::Mating_Plano(molecule* mom, molecule* dad, TRandom3* gRandom){

  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();
  double pos_CM = 0;

  int flag = 0;
  int nr_atoms = 0;

  //gRandom = new TRandom3(0);
  int dir = (int)gRandom -> Uniform(0,3);

  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++)
    positions[i][j] = 0;
  }

  
  //posicao centro de massa 
  for(int i = 0; i < N_atoms; i++)
    pos_CM += (pos_mom[i][dir] + pos_dad[i][dir])/2;
  pos_CM = pos_CM/N_atoms; 


  //choose every mom atom above CM
  for (int i = 0; i < N_atoms; i++){


    if(pos_mom[i][dir] > pos_CM){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_mom[i][j];
      nr_atoms ++;
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
        }
      }
      if(nr_atoms== N_atoms)
        break;
    }
  }

}


