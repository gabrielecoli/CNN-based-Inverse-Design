#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mc.h"

int main(int argc, char **argv) {
  vector *par;
  box B;
  unsigned long long int Nattempts_dis = 0, Naccepted_dis = 0, n;
  unsigned long long int Nattempts_vol = 0, Naccepted_vol = 0;
  #ifdef NPT
  double density;
  FILE *den;
  #endif
  double acceptance_dis, acceptance_vol, Dr, energy;
  double av_energy = 0.0, av_density = 0.0;
  int Nmeasures = 0;
  int i;
  #ifdef NVT
  double g[Nbins] = {0.0};
  int Ng = 0;
  #endif
  double DeltaPos, DeltaL, minL;
  vector scaling;
  FILE *en, *acc, *info, *fpsiminfo, *fpfitness, *fpflag;
  char file_string[150];
  // FOR BOPS
  int m;
  int *Nb;
  neighbor **neighbors;
  double *radius;
  double **chi;
  double av_neighbors;
  double Chi[m_MAX+1];
  FILE *op;
  // END FOR BOPS

  scaling.x = 1.0;
  scaling.y = 1.0;
  scaling.z = 1.0;


  sprintf(file_string, "siminfo.dat");
  if((fpsiminfo = fopen(file_string,"r"))==NULL){
      printf("ERROR: File with simulation parameters could not be opened.\n");
      exit(0);
  }
  fscanf(fpsiminfo, "%lf %lf %lf", &T, &InitDensity, &sigma1);
  init_param(&B, &DeltaPos, &DeltaL, &Dr);

  par = (vector *) malloc(Npart * sizeof(vector));
  if(par == NULL) {
    printf("\nERROR: Malloc of par failed.\n");
  }

  // FOR BOPS
  Nb = (int *) malloc(Npart * sizeof(int));
  neighbors = (neighbor **) malloc(Npart * sizeof(neighbor*));
  for(i=0; i<Npart; i++) {
    neighbors[i] = (neighbor *) malloc((MAX_NEIGH) * sizeof(neighbor));
  }
  radius = (double *) malloc(Npart * sizeof(double));
  chi = (double **) malloc(Npart * sizeof(double*));
  for(i=0; i<Npart; i++) {
    chi[i] = (double *) malloc((m_MAX+1) * sizeof(double));
  }
  // END FOR BOPS



  if(InitConf==0) {
    init_rand_conf(par, B);
  }
  else if(InitConf==1) {
    init_sc_conf(par, B);
  }
  else if(InitConf==2){
    init_hex_conf(par, B);
  }
  else {
    read_init_conf(par, &B);
  }
  n = 0;
  print_snapshot(par, B, n, argv);

  sprintf(file_string, "info.txt");
  if((info = fopen(file_string,"w"))==NULL) {
    printf("\nError while opening file [info.txt]\n");
    exit(0);
  }

  // print parameters of the simulation
  #ifdef NVT
  printf("\n**********        MC NVT-simulation        *********\n");
  fprintf(info, "\n**********        MC NVT-simulation        *********\n");
  #endif
  #ifdef NPT
  printf("\n**********        MC NPT-simulation        *********\n");
  fprintf(info, "\n**********        MC NPT-simulation        *********\n");
  #endif
  printf("\n**********  Parameters of the simulation  **********\n");
  fprintf(info, "\n**********  Parameters of the simulation  **********\n");

  printf("Simulation seed =                           %d\n", seed);
  fprintf(info, "Simulation seed =                           %d\n", seed);
  printf("Number of particles =                       %d\n", Npart);
  fprintf(info, "Number of particles =                       %d\n", Npart);
  printf("T =                                         %lg\n", T);
  fprintf(info, "T =                                         %lg\n", T);
  printf("BetaP =                                     %lg\n", betaP);
  fprintf(info, "BetaP =                                     %lg\n", betaP);
  printf("Number of MC equilibration iterations =     %llu\n", Nequilibrations);
  fprintf(info, "Number of MC equilibration iterations =     %llu\n", Nequilibrations);
  printf("Number of MC iterations in equilibrium =    %llu\n", Niterations);
  fprintf(info, "Number of MC iterations in equilibrium =    %llu\n", Niterations);
  printf("Density =                                   %lg\n", InitDensity);
  fprintf(info, "Density =                                   %lg\n", InitDensity);
  printf("Lx =                                        %lg\n", B.x);
  fprintf(info, "Lx =                                        %lg\n", B.x);
  printf("Ly =                                        %lg\n", B.y);
  fprintf(info, "Ly =                                        %lg\n", B.y);
  printf("Lz =                                        %lg\n", B.z);
  fprintf(info, "Lz =                                        %lg\n", B.z);
  printf("Maximum displacement =                      %lg\n", DeltaPos);
  fprintf(info, "Maximum displacement =                      %lg\n", DeltaPos);
  printf("Maximum volume change =                     %lg\n", DeltaL);
  fprintf(info, "Maximum volume change =                     %lg\n", DeltaL);
  printf("sigma1 =                                    %lg\n", sigma1);
  fprintf(info, "sigma1 =                                    %lg\n", sigma1);


  // Equilibration
  sprintf(file_string, "energy.dat");
  if((en = fopen(file_string,"w"))==NULL) {
    printf("\nError while opening file [energy.dat]\n");
    exit(0);
  }
  sprintf(file_string, "boop.dat");
  if((op = fopen(file_string,"w"))==NULL) {
    printf("\nError while opening file [boop.dat]\n");
    exit(0);
  }
  #ifdef NPT
  sprintf(file_string, "density.dat");
  if((den = fopen(file_string,"w"))==NULL) {
    printf("\nError while opening file [density.dat]\n");
    exit(0);
  }
  #endif
  sprintf(file_string, "acceptance.dat");
  if((acc = fopen(file_string,"w"))==NULL) {
    printf("\nError while opening file [acceptance.dat]\n");
    exit(0);
  }
  for(n=0; n<Nequilibrations; n++) {
    for(i=0; i<Npart; i++){
      displacement_move(par, B, DeltaPos, &Nattempts_dis, &Naccepted_dis);
    }
    #ifdef NPT
    volume_move(par, &B, DeltaL, &Nattempts_vol, &Naccepted_vol);
    #endif
    if(n%Nadjust==0) {
      // compute acceptance
      acceptance_dis = (double) Naccepted_dis/Nattempts_dis;
      acceptance_vol = (double) Naccepted_vol/Nattempts_vol;
      fprintf(acc, "%llu %lf %lf %lf %lf\n", n, DeltaPos, acceptance_dis, DeltaL, acceptance_vol);
      fflush(acc);
      // adjust the moves displacements
      if(acceptance_dis>(TARGET_ACC_DIS+TOL)) {
        DeltaPos *= 1.1;
      }
      else if(acceptance_dis<(TARGET_ACC_DIS-TOL)){
        DeltaPos /= 1.1;
      }
      if(DeltaPos > B.x/4) {
        DeltaPos = B.x/4;
      }
      if(acceptance_vol>(TARGET_ACC_VOL+TOL)) {
        DeltaL *= 1.05;
      }
      else if(acceptance_vol<(TARGET_ACC_VOL-TOL)){
        DeltaL /= 1.05;
      }
      minL = B.x;
      if(B.y < minL) {
        minL = B.y;
      }
      if(sqrt(DeltaL) > minL/3) {
        DeltaL = minL*minL/9;
      }
      // reset scceptances to zero
      Nattempts_dis = 0;
      Naccepted_dis = 0;
      Nattempts_vol = 0;
      Naccepted_vol = 0;
    }
    if(n%Nprint==0) {
      #ifdef NPT
      density = Npart/B.V;
      fprintf(den, "%llu %lf\n", n, density);
      fflush(den);
      #endif
      energy = energy_system(par, B, scaling)/Npart;
      fprintf(en, "%llu %lf\n", n, energy);
      fflush(en);
      //print_snapshot(par, B, n, argv);
      // FOR BOPS
      for(i=0; i<Npart; i++) {
        neighbors_sann(par, B, Nb, radius, neighbors, chi, i);
        bops(par, Nb, chi, neighbors, i);
      }
      av_neighbors = 0.0;
      for(m=1; m<m_MAX+1; m++) {
        Chi[m] = 0.0;
      }
      for(i=0; i<Npart; i++) {
        av_neighbors += (double)Nb[i]/Npart;
        for(m=1; m<m_MAX+1; m++) {
          Chi[m] += chi[i][m]/Npart;
        }
      }
      fprintf(op, "%llu %lf ", n, av_neighbors);
      for(m=1; m<m_MAX+1; m++) {
        fprintf(op, "%lf ", Chi[m]);
      }
      fprintf(op, "\n");
      fflush(op);
      // END FOR BOPS
    }
  }
  fclose(acc);
  // reset scceptances to zero
  Nattempts_dis = 0;
  Naccepted_dis = 0;
  Nattempts_vol = 0;
  Naccepted_vol = 0;

  // production
  for(n=1; n<=Niterations; n++) {
    for(i=0; i<Npart; i++){
      displacement_move(par, B, DeltaPos, &Nattempts_dis, &Naccepted_dis);
    }
    #ifdef NPT
    volume_move(par, &B, DeltaL, &Nattempts_vol, &Naccepted_vol);
    #endif
    if(n%Nprint==0) {
      Nmeasures += 1;
      #ifdef NVT
      gr(g, Dr, par, B, &Ng);
      av_density += InitDensity;
      #endif
      #ifdef NPT
      density = Npart/B.V;
      av_density += density;
      fprintf(den, "%llu %lf\n", n+Nequilibrations, density);
      fflush(den);
      #endif
      energy = energy_system(par, B, scaling)/Npart;
      av_energy += energy;
      fprintf(en, "%llu %lf\n", n+Nequilibrations, energy);
      fflush(en);
      print_snapshot(par, B, n, argv);
      // FOR BOPS
      for(i=0; i<Npart; i++) {
        neighbors_sann(par, B, Nb, radius, neighbors, chi, i);
        bops(par, Nb, chi, neighbors, i);
      }
      av_neighbors = 0.0;
      for(m=1; m<m_MAX+1; m++) {
        Chi[m] = 0.0;
      }
      for(i=0; i<Npart; i++) {
        av_neighbors += (double)Nb[i]/Npart;
        for(m=1; m<m_MAX+1; m++) {
          Chi[m] += chi[i][m]/Npart;
        }
      }
      fprintf(op, "%llu %lf ", n+Nequilibrations, av_neighbors);
      for(m=1; m<m_MAX+1; m++) {
        fprintf(op, "%lf ", Chi[m]);
      }
      fprintf(op, "\n");
      fflush(op);
      // END FOR BOPS
    }
  }
  fclose(en);
  #ifdef NPT
  fclose(den);
  #endif
  #ifdef NVT
  // normalize_gr(g, Dr, Ng, InitDensity, argv);
  #endif
  // normalize averages
  av_density /= Nmeasures;
  av_energy /= Nmeasures;
  // compute acceptance
  acceptance_dis = (double) Naccepted_dis/Nattempts_dis;
  acceptance_vol = (double) Naccepted_vol/Nattempts_vol;
  // print reslts on screen
  printf("Acceptance displacement moves   = %lg\n", acceptance_dis);
  fprintf(info, "Acceptance displacement moves   = %lg\n", acceptance_dis);
  printf("Acceptance volume moves         = %lg\n", acceptance_vol);
  fprintf(info, "Acceptance volume moves         = %lg\n", acceptance_vol);
  printf("Average energy per particle     = %lg\n", av_energy);
  fprintf(info, "Average energy per particle     = %lg\n", av_energy);
  printf("Average density                 = %lg\n", av_density);
  fprintf(info, "Average density                 = %lg\n", av_density);
  fclose(info);

  /*
  sprintf(file_string, "../generation%s/sample%s/fitness.dat", argv[1], argv[2]);
  if((fpfitness = fopen(file_string,"w+"))==NULL){
    printf("ERROR: File with fitness value could not be opened.\n");
    exit(0);
  }
  fprintf(fpfitness, "%lf\n", 100*SQR(Chi[SymmetryIndex] - Target_QC));
  fflush(fpfitness);
  fclose(fpfitness);
  */

  free(par);
  // FOR BOPS
  free(Nb);
  free(radius);
  for(i=0; i<Npart; i++) {
    free(neighbors[i]);
    free(chi[i]);
  }
  free(neighbors);
  free(chi);
  fclose(op);
  // END FOR BOPS

  
  sprintf(file_string, "flag.dat");
  if((fpflag = fopen(file_string,"w"))==NULL){
    printf("ERROR: Flag did not work.\n");
    exit(0);
  }
  fprintf(fpflag, "1\n");
  fflush(fpflag);
  fclose(fpflag);
  

  return 0;
}


void init_param(box *b, double *DeltaPos, double *DeltaL, double *Dr) {
  double V, L, tempdensity;
  FILE *f;
  char buf[200];
  char temp[200];
  char *res;

  // open input file
  f = fopen("in.dat", "r");
  if(f==NULL) {
    printf("\nERROR: file [in.dat] not found!\n");
    exit(0);
  }
  // read input file
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %d", temp, &Npart);
  if(strcmp(temp, "Npart:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [Npart:]\n", temp);
    exit(1);
  }
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %lf", temp, &tempdensity);
  if(strcmp(temp, "InitDensity:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [InitDensity:]\n", temp);
    exit(1);
  }
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %llu", temp, &Nequilibrations);
  if(strcmp(temp, "Nequilibrations:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [Nequilibrations:]\n", temp);
    exit(1);
  }
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %llu", temp, &Niterations);
  if(strcmp(temp, "Niterations:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [Niterations:]\n", temp);
    exit(1);
  }
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %d", temp, &Nprint);
  if(strcmp(temp, "Nprint:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [Nprint:]\n", temp);
    exit(1);
  }
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %d", temp, &seed);
  if(strcmp(temp, "seed:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [seed:]\n", temp);
    exit(1);
  }
  // line
  res = fgets(buf, 200, f);
  sscanf(buf, "%s %d", temp, &InitConf);
  if(strcmp(temp, "InitConf:")!=0) {
    printf("ERROR: while reading input file \nfound [%s] instead of [InitConf:]\n", temp);
    exit(1);
  }
  if(InitConf > 2 && InitConf<0) {
    printf("ERROR: non-valid value for InitConf(%d)!\nMust be 0, 1 or 2.\n", InitConf);
    exit(1);
  }
  fclose(f);

  // init random number generator
  if(seed == 0) {
    seed = time(0);
  }
  srand48(seed);
  // initialize cubic box
  V = Npart/InitDensity;
  if(InitConf ==2) {
    L = sqrt(2.0*V/sqrt(3.0));
    if(rc>L*0.5) {
      printf("\nERROR: rc larger than half of the box lenght [L = %lg]!\n", L);
      printf("\n       Increase Npart or reduce density.\n");
      exit(0);
    }
    b->x = L;
    b->y = L*sqrt(3.0)*0.5;
    b->z = 1.0;
    b->V = V;
  }
  else {
    L = sqrt(V);
    b->x = L;
    b->y = L;
    b->z = 1.0;
    b->V = V;
  }
  *Dr = 0.5*b->y/Nbins;
  *DeltaPos = 0.1;
  *DeltaL = V/100.0;
  beta = 1.0/T;
}


void init_rand_conf(vector *par, box b) {
  int i;
  vector scale;

  scale = setVector(1, 1, 1);
  for(i=0; i<Npart; i++) {
    // random position
    par[i] = setVector(ran()*b.x,ran()*b.y,0.0);
  }
  printf("\nInitial random configuration with %d particles succesfully generated!\n", Npart);
}

void init_sc_conf(vector *par, box b) {
  int i, j, count, Nlattice;
  double x, y, z, spacing, diff;

  Nlattice = sqrt(Npart);
  diff = sqrt(Npart) - Nlattice;
  if(diff>0.0000001) {
    Nlattice += 1;
  }
  spacing = b.x/Nlattice;
  count = 0;

  x = 0.0;
  for(i=0; i<Nlattice; i++) {
    y = 0.0;
    for(j=0; j<Nlattice; j++) {
      z = 0.0;
	    if(count<Npart) {
	      par[count] = setVector(x,y,z);
	      count++;
      }
      y += spacing;
    }
    x += spacing;
  }
  printf("\nInitialization on a sc lattice: %d particles placed\n", count);
}

void init_hex_conf(vector *par, box b) {
  int i, j, count, Nlattice;
  double x, y, z, spacing, sigma, diff;
  box tempb;
  vector scale;

  Nlattice = sqrt(Npart);
  diff = sqrt(Npart) - Nlattice;
  if(diff>0.0000001) {
    Nlattice += 1;
  }
  spacing = 1.0;
  sigma = 1.0;
  count = 0;

  for(i=0; i<Nlattice; i++) {
    for(j=0; j<Nlattice; j++) {
	    if(count<Npart) {
        if(j%2 == 0) {
          x = i*spacing*sigma;
        }
        else {
          x = (i+0.5)*spacing*sigma;
        }
        y = j*spacing*sigma*sqrt(3.0)*0.5;
        z = 0.0;
	      par[count] = setVector(x,y,z);
	      count++;
      }
    }
  }
  tempb.x = Nlattice*spacing*sigma;
  tempb.y = Nlattice*spacing*sigma*sqrt(3.0)*0.5;
  tempb.z = 1.0;
  tempb.V = tempb.x*tempb.y*tempb.z;
  // rescale particle positions
  scale = setVector(b.x/tempb.x, b.y/tempb.y, 1);
  for(i=0; i<Npart; i++) {
    par[i].x *= scale.x;
    par[i].y *= scale.y;
  }
  // rescale box
  tempb.x *= scale.x;
  tempb.y *= scale.y;
  tempb.V = tempb.x*tempb.y*tempb.z;

  printf("\nInitialization on a hexagonal lattice: %d particles placed\n", count);
}

void read_init_conf(vector *par, box *b) {
  int i, N;
  double temp;
  box temp_box;
  FILE *in;

  if((in = fopen("init_conf.dat","r"))==NULL) {
    printf("Error while opening file\n");
    exit(0);
  }
  fscanf(in,"%d\n", &N);
  fscanf(in,"%lf %lf\n", &temp, &temp_box.x);
  fscanf(in,"%lf %lf\n", &temp, &temp_box.y);
  fscanf(in,"%lf %lf\n", &temp, &temp_box.z);
  b->x = temp_box.x;
  b->y = temp_box.y;
  b->z = temp_box.z;
  for(i=0; i<Npart; i++) {
    fscanf(in,"%lf %lf %lf %lf\n", &par[i].x, &par[i].y, &par[i].z, &temp);
  }
  fclose(in);

  printf("\nInitial configuration with %d particles succesfully loaded from file!\n", Npart);
}

void print_snapshot(vector *par, box b, unsigned long long int n, char ** argv) {
  int i;
  char file_name[80];
  FILE *out;

  sprintf(file_name, "conf/conf_%06llu.dat", n);
  if((out = fopen(file_name,"w"))==NULL) {
    printf("Error while opening file\n");
    exit(0);
  }
  fprintf(out,"%d\n", Npart);
  //fprintf(out,"%lg %lg %lg\n", b.x, b.y, b.z);
  fprintf(out,"0 0 0\n");
  fprintf(out,"%lf 0 0\n", b.x);
  fprintf(out,"0 %lf 0\n", b.y);
  fprintf(out,"0 0 %lf\n", b.z);
  for(i=0; i<Npart; i++) {
    fprintf(out,"%lf %lf %lf 1.00 1\n", par[i].x, par[i].y, par[i].z);
  }
  fclose(out);
}


void volume_move(vector *par, box *b, double DeltaL, unsigned long long int *Nattempts, unsigned long long int *Naccepted) {
  double gamma, ratio, newV, newVolumaccio;
  vector scalen, scaleo;
  int i;
  double eno, enn;
  double r, p;


  *Nattempts += 1;
  // random change of volume
  newVolumaccio = b->V + DeltaL*(ran()-0.5);
  ratio = b->x/b->y;
  //printf("ratio is: %lf\n", ratio);
  gamma = newVolumaccio/b->V;
  scalen = setVector(pow(gamma,0.5), pow(gamma,0.5), 1);
  // compute new volume
  newV = scalen.x*b->x*scalen.y*b->y*scalen.z*b->z;
  if(fabs(newV-newVolumaccio)>0.001){
    printf("ERROR! VOLUMACCIO PROBLEM");
  }
  // energy of old configuration
  scaleo = setVector(1, 1, 1);
  eno = 0.0;
  for(i=0; i<Npart; i++) {
    //printf("newV = %lf, oldV = %lf, calculating old energy...\n", newV, b->V);
    eno += energy_particle(par[i], i, par, *b, scaleo, i+1);
  }
  // energy of new configuration
  enn = 0.0;
  for(i=0; i<Npart; i++) {
    //printf("newV = %lf, oldV = %lf, calculating new energy...\n", newV, b->V);
    enn += energy_particle(par[i], i, par, *b, scalen, i+1);
  }
  r = ran();
  p = exp(-beta*(enn-eno)-betaP*(newV-b->V))*pow(newV/b->V, Npart);
  if(r<p) {
    // accept
    *Naccepted += 1;
    // update positions and volume
    rescale(par, b, scalen, newV);
  }
}


void displacement_move(vector *par, box b, double DeltaPos, unsigned long long int *Nattempts, unsigned long long int *Naccepted) {
  vector dr, newpos, scale;
  double eno, enn, r, p;
  int i;

  //printf("Displacement move!!!");
  *Nattempts += 1;
  // pick a particle at random
  i = lrand48() % Npart;
  // random displacement
  dr = setVector((ran()-0.5)*DeltaPos,(ran()-0.5)*DeltaPos,0.0);
  newpos = sumVectors(par[i], dr);
  // setting new volume equal to old volume
  scale = setVector(1, 1, 1);
  // energy of new configuration
  enn = energy_particle(newpos, i, par, b, scale, 0);
  // energy of old configuration
  eno = energy_particle(par[i], i, par, b, scale, 0);
  //metropolis check
  r = ran();
  p = exp(-beta*(enn-eno));
  if(r<p) {
    //accept
    *Naccepted += 1;
    // put particle in simulation box
    if(newpos.x < 0.0)
      newpos.x += b.x;
    else if(newpos.x >= b.x)
      newpos.x-=b.x;

    if(newpos.y < 0.0)
      newpos.y += b.y;
    else if(newpos.y >= b.y)
      newpos.y -= b.y;

    if(newpos.z < 0.0)
      newpos.z += b.z;
    else if(newpos.z >= b.z)
      newpos.z -= b.z;
    //update new position
    par[i] = setVector(newpos.x, newpos.y, newpos.z);
  }
}

void rescale(vector *par, box *b, vector scaling, double newV) {
  int i;

  // rescale particle positions
  for(i=0; i<Npart; i++) {
    par[i].x *= scaling.x;
    par[i].y *= scaling.y;
    par[i].z *= scaling.z;
  }
  // rescale box
  b->V = newV;
  b->x *= scaling.x;
  b->y *= scaling.y;
  b->z *= scaling.z;
}



double energy_particle(vector newpos, int i, vector *par, box b, vector scaling, int k) {
  int j;
  vector d;
  double r;
  double energy;

  energy = 0.0;
  for(j=k; j<Npart; j++) {
    if(j != i) {
      // distance
      d.x = scaling.x*newpos.x - scaling.x*par[j].x;
      d.y = scaling.y*newpos.y - scaling.y*par[j].y;
      d.z = scaling.z*newpos.z - scaling.z*par[j].z;
      // periodic boundary conditions
      d = pbc(d, b, scaling);
      // radial distance
      r = sqrt(dotProduct(d,d));
      // add neighbour
      if(r < b.x*0.5) {
        // 2-body energy contribution
        energy += epsilon*(pow(1./r, 14) + 0.5*(1.-tanh(kappa*(r-sigma1))));
      }
    }
  }
  return energy;
}

int overlap_check(vector newpos, int i, vector *par, box b, vector scaling) {
  int j;
  vector d;
  double r2;

  for(j=0; j<Npart; j++) {
    if(j != i) {
      // distance
      d.x = scaling.x*newpos.x - scaling.x*par[j].x;
      d.y = scaling.y*newpos.y - scaling.y*par[j].y;
      d.z = scaling.z*newpos.z - scaling.z*par[j].z;
      // periodic boundary conditions
      d = pbc(d, b, scaling);

      r2 = dotProduct(d,d);
      //check
      if(r2<=rc*rc) {
        // overlap
        return 1;
      }
    }
  }
  return 0;
}

int partial_overlap_check(vector newpos, int i, vector *par, box b, vector scaling) {
  int j;
  vector d;
  double r2;

  for(j=0; j<i; j++) {
    // distance
    d.x = scaling.x*newpos.x - scaling.x*par[j].x;
    d.y = scaling.y*newpos.y - scaling.y*par[j].y;
    d.z = scaling.z*newpos.z - scaling.z*par[j].z;
    // periodic boundary conditions
    d = pbc(d, b, scaling);

    r2 = dotProduct(d,d);
    //check
    if(r2<=rc*rc) {
      // overlap
      return 1;
    }
  }
  return 0;
}


// sampling

double energy_system(vector *par, box b, vector scaling) {
  int i;
  double energy;

  energy = 0.0;
  for(i=0; i<Npart; i++){
    energy += energy_particle(par[i], i, par, b, scaling, i+1);
  }
  return energy;
}


void gr(double *g, double Dr, vector *par, box b, int *N) {
  int i, j, k;
  vector d, scaling;
  double r;

  scaling = setVector(1,1,1);
  *N += 1;
  for(i=0; i<Npart-1; i++) {
    for(j=i+1; j<Npart; j++) {
      // distance
      d.x = par[i].x - par[j].x;
      d.y = par[i].y - par[j].y;
      d.z = par[i].z - par[j].z;
      // periodic boundary conditions
      d = pbc(d, b, scaling);
      // radial distance
      r = sqrt(dotProduct(d,d));
      // only within half of the box lenght
      if(r <= b.y*0.5) {
        k = r/Dr;
      	g[k] += 2.0;
      }
    }
  }
}

void normalize_gr(double *g, double Dr, int Ng, double rho, char ** argv) {
  int i;
  double r, vb, Nid;
  char file_name[80];
  FILE *out;

  sprintf(file_name, "g_rho%.3f.dat", InitDensity);
  if((out = fopen(file_name, "w"))==NULL) {
    printf("Error while opening file\n");
    exit(0);
  }
  for(i=0; i<Nbins; i++) {
    // distance r
    r = Dr*(i+0.5);
    // volume between bin i+1 and i
    vb = (pow(i+1,2)-pow(i,2))*pow(Dr,2);
    // number of ideal gas particle in vb
    Nid = M_PI*vb*rho;
    // normalization of g
    g[i] = g[i]/(Ng*Npart*Nid);
    fprintf(out, "%lf %lf\n", r, g[i]);
  }
  fclose(out);
}

// vectors

vector pbc(vector dist, box b, vector scaling) {
  vector d;
  // apply periodic boundary conditions
  d.x = dist.x - b.x*scaling.x*rint(dist.x/(b.x*scaling.x));
  d.y = dist.y - b.y*scaling.y*rint(dist.y/(b.y*scaling.y));
  d.z = dist.z - b.z*scaling.z*rint(dist.z/(b.z*scaling.z));

  return d;
}

double ran() {
  double r;
  // random number in (0,1)
  r = (double) lrand48()/RAND_MAX;
  return r;
}

vector setVector(double x, double y, double z) {
  vector result;
  result.x = x;
  result.y = y;
  result.z = z;
  return result;
}

vector sumVectors(vector a, vector b) {
  vector result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  result.z = a.z + b.z;
  return result;
}

vector distanceVector(vector a, vector b) {
  vector d;
  // distance a-b
  d.x = a.x-b.x;
  d.y = a.y-b.y;
  d.z = a.z-b.z;
  return d;
}

double dotProduct(vector a, vector b) {
  double dot;
  dot = a.x*b.x + a.y*b.y + a.z*b.z;
  return dot;
}


int compare(const void *nb1, const void *nb2) {
  neighbor *pnb1 = (neighbor*) nb1;
  neighbor *pnb2 = (neighbor*) nb2;
  // if r1 < r2 return -1
  if(pnb1->r < pnb2->r) return -1;
  return 1;
}

void neighbors_sann(vector *par, box b, int *Nb, double *radius, neighbor **neighbors, double **chi, int i) {
  int j, m, k;
  double dx, dy, dz, dr;
  double sum;

  Nb[i] = 0;
  for(j=0; j<Npart; j++) {
    if(j != i) {
      // distance
      dx = par[j].x - par[i].x;
      dy = par[j].y - par[i].y;
      dz = par[j].z - par[i].z;
      // boundary conditions
      dx -= b.x*rint(dx/b.x);
      dy -= b.y*rint(dy/b.y);
      dz -= b.z*rint(dz/b.z);
      // distance
      dr = sqrt(dx*dx+dy*dy+dz*dz);
      // check
      if(dr<=rc_sann) {
        // update number of neighbors
        Nb[i] = Nb[i]+1;
        if(Nb[i] > MAX_NEIGH) {
          printf("ERROR: increase the maximum number of neighbors (MAX_NEIGH)!\n");
          exit(0);
        }
        // add j as a neighbor of i
        neighbors[i][Nb[i]-1].id = j;
        // add distance between i and j to the list
        neighbors[i][Nb[i]-1].r = dr;
        neighbors[i][Nb[i]-1].x = dx;
        neighbors[i][Nb[i]-1].y = dy;
        neighbors[i][Nb[i]-1].z = dz;
      }
    }
  }
  // Step 2: sort neighbors of particle iaccording to rheir distance
  // if there were not enough neighbors report an error
  if(Nb[i] < 3) {
    printf("ERROR: only %d neighbours found for particle %d within given cutoff!\n", Nb[i], i);
    exit(0);
  }
  qsort(neighbors[i], Nb[i], sizeof(neighbor), compare);
  // Step 3: start with m = 3
  m = 3;
  // compute R(m) as in Eq.3
  sum = 0.0;
  for(k=0; k<m; k++) {
    sum = sum + neighbors[i][k].r;
  }
  radius[i] = sum/(k-2.0);
  // Step 4/5: include further neighbors until finished
  //(SANN radius smaller than next neighbor distance)
  while((k<Nb[i]) && radius[i]>neighbors[i][k].r) {
    sum += neighbors[i][k].r;
    radius[i] = sum/(k-1.0);
    k++;
  }
  // if not enough neighbors to converge report an error
  if(k == Nb[i]) {
      printf("ERROR: not enough (%d) neighbours for particle %d within given cutoff.\nSANN could not converge!\n", Nb[i], i);
      exit(0);
    }
  // Step 6: save number of neighbors and selected neighbors
  Nb[i] = k;
  radius[i] = radius[i];
}

void bops(vector *par, int *Nb, double **chi, neighbor **neighbors, int i) {
  int j, m;
  double theta;
  complex temp[m_MAX+1];

  for(m=1; m<m_MAX+1; m++) {
    temp[m].Re = 0.0;
    temp[m].Im = 0.0;
  }

  if(Nb[i]>0) {
    for(j=0; j<Nb[i]; j++) {
      theta = atan2(neighbors[i][j].y, neighbors[i][j].x);
      for(m=1; m<m_MAX+1; m++) {
        temp[m].Re += cos(m*theta);
        temp[m].Im += sin(m*theta);
      }
    }
    for(m=1; m<m_MAX+1; m++) {
      temp[m].Re /= Nb[i];
      temp[m].Im /= Nb[i];
      chi[i][m] = temp[m].Re*temp[m].Re + temp[m].Im*temp[m].Im;
    }
  }
  else {
    for(m=1; m<m_MAX+1; m++) {
      chi[i][m] = 0.0;
    }
  }
}
