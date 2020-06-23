#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#define DT              1.0e-2
//#define OUTDT           1.0e-0 /*changed 0130*/
#define TMAX            1.0e+2///////////////
#define MAXITER         ((int)((TMAX)/(DT)))
#define OUTGAP          ((int)((MAXITER)/1000))

#define D               1.0e-2  //um^2/s
#ifndef RAMDA 
#define RAMDA           0.0     //Sub-diffusion Parameter in Ornstein-Uhlenbeck process 
#endif
#define Lx              1.0e0   //initial length
#define Ly              1.0e0   //initial length
#define lx              1.0e0   //initial length
#define ly              1.0e0   //initial length
#define l0              Lx/2.0
#define T0              (l0*l0)/(9.0*2.0*D)
#define K               0.000  //Parameter : Strain/Diffusion 
#define Epx             K/(4*T0)  //*100%/s 
#define Epy             Epx //0.0e0   //*100%/s
#define si              sqrt(2.0*D*DT)
#define SIGMA           sqrt(2.0*D)
#define NUM_PARTICLE    (int)(10000)
#define THETAX          0.0
#define THETAY          0.0
#define RAMDAX          RAMDA 
#define RAMDAY          RAMDA 
#define ALPHA           1.0
#define H               ALPHA/2.0
double randRange(double min, double max){return ((rand()+1.0)/((RAND_MAX+1.0)/(max - min))+min);}
double rand_normal( double mu, double sigma){
    double z=sqrt(-2.0*log(randRange(0.0,1.0)))*sin(2.0*M_PI*randRange(0.0,1.0)); 
    return (mu + sigma*z);
}
double Strain(double L0, double de, double t){
  return(L0 * (1.0 + de*t) );
}
class Particle
{
  public:
    Particle();
    ~Particle();

  public:
    void    Update(int iter);
    void    Output(int iter);
  private:
    double  *x_, *v_, *xi_, *xgl_;
    double  *L_;
    int     numParticle_;
    
};

Particle::Particle()
{
  numParticle_ = NUM_PARTICLE; 

  x_ = new double[2*numParticle_];
  v_ = new double[2*numParticle_];
  xi_ = new double[2*numParticle_];
  xgl_ = new double[2*numParticle_];
  L_ = new double[2];
  srand(1);
  for(int idParticle = 0; idParticle < numParticle_; idParticle++){
    x_[2*idParticle + 0] = 0.00000000;
    x_[2*idParticle + 1] = 0.00000000;
    xi_[2*idParticle + 0]= x_[2*idParticle + 0];//Initial position
    xi_[2*idParticle + 1]= x_[2*idParticle + 1];//Initial position   
    xgl_[2*idParticle + 0] = x_[2*idParticle + 0];
    xgl_[2*idParticle + 1] = x_[2*idParticle + 1];
    L_[0] = Lx; 
    L_[1] = Ly;
  }
}

Particle::~Particle()
{
  delete [] x_;
  delete [] v_;
  delete [] xi_;
  delete [] xgl_;
  delete [] L_;
}

void Particle::Update(int iter)
{
 // double sum=0;
 numParticle_ = NUM_PARTICLE; 
  double r[2]; 
  double e=1.0/2.0*( pow((double) iter*DT,(2.0*H)) + pow(DT,(2.0*H)) - pow((double) DT*(iter-1.0),(2.0*H)) );
//printf("H=%.4f E=%.4f\n",H,e);
//  double delx[2];
  double delxgl[2];
  for(int idParticle = 0; idParticle < numParticle_; idParticle++){
    r[0]=rand_normal(0.000,sqrt(e));
    r[1]=rand_normal(0.000,sqrt(e));
//    delx[0]    = -RAMDAX*(x_[2*idParticle +0] - THETAX)*DT + SIGMA*r[0];
//    delx[1]    = -RAMDAY*(x_[2*idParticle +1] - THETAY)*DT + SIGMA*r[1];
    delxgl[0]  = -RAMDAX*(xgl_[2*idParticle +0] - THETAX)*DT + SIGMA*r[0];
    delxgl[1]  = -RAMDAY*(xgl_[2*idParticle +1] - THETAY)*DT + SIGMA*r[1];
    x_[2*idParticle + 0] += delxgl[0];
    x_[2*idParticle + 1] += delxgl[1];
    xgl_[2*idParticle + 0] += delxgl[0];
    xgl_[2*idParticle + 1] += delxgl[1];
    if (x_[2*idParticle + 0] >= L_[0]/2.0){
            x_[2*idParticle + 0] = x_[2*idParticle + 0] - L_[0]; 
          
     }else if (x_[2*idParticle + 0] <= - L_[0]/2.0){
            x_[2*idParticle + 0] = x_[2*idParticle + 0] + L_[0]; 
    }

    if (x_[2*idParticle + 1] >= L_[1]/2.0){
            x_[2*idParticle + 1] = x_[2*idParticle + 1] - L_[1]; 
     }
    else if (x_[2*idParticle + 1] <= - L_[1]/2.0){
            x_[2*idParticle + 1] = x_[2*idParticle + 1] + L_[1]; 
    }
 }
}
void Particle::Output(int iter)
{
/*calculation for msd*/
  double dx,dy;
  double dr2;
  double sum2 = 0.0;
  double A1   = 0.0;
  double A2   = 0.0;
  double A3   = 0.0;
//  double num  = 0.0;
  double msd  = 0.0;
  double msd4 = 0.0;
  if(iter%OUTGAP != 0) return;
  for (int idParticle = 0; idParticle < numParticle_; idParticle++){
    dx  = xgl_[2*idParticle + 0]-xi_[2*idParticle + 0];
    dy  = xgl_[2*idParticle + 1]-xi_[2*idParticle + 1];
   //dx  = x_[2*idParticle + 0]-xi_[2*idParticle + 0];
   //dy  = x_[2*idParticle + 1]-xi_[2*idParticle + 1];
    dr2 = dx*dx + dy*dy;
    sum2 += dr2;
    A1 += (dx*dx*dx*dx + dy*dy*dy*dy);
    A2 += (dx*dx);
    A3 += (dy*dy);
/*/  if ( (xgl_[2*idParticle + 1]<L_[1]/2.0 && xgl_[2*idParticle + 1]>-L_[1]/2.0)\
      && (xgl_[2*idParticle + 0]<L_[0]/2.0 && xgl_[2*idParticle + 0]>-L_[0]/2.0)){
      num+=1;
    }*/
  }
  msd   = sum2 / (double) numParticle_;                                                                 //MSD
  msd4  = A1 / (double) numParticle_ + 2.0* (A2 / (double) numParticle_) *(A3 / (double) numParticle_) ;//forth moment of dr ensunbled with N
  //double alpha = 3.0*msd4 / (5.0*msd*msd) - 1.0; 
  double alpha = msd4 / (2.0*msd*msd) - 1.0; 
///////////////////////////////////
  FILE *fp;
  char filename[256];//for trajectories
  FILE *gp;
  char filename2[256];//for displaying the wall
  FILE *hp;
  char wallfile[256];//for MSDs
//////////////////////////////////////////////
  sprintf(filename, "result/result%05d.dat",iter/OUTGAP);
  if((fp = fopen(filename, "w")) == NULL){printf("FAILED TO OPEN 1st FILE.\n"); exit(1);};
    for(int idParticle = 0; idParticle < numParticle_; idParticle++){
      fprintf(fp, "%15e %15e %15e %15e\n", x_[2*idParticle + 0], x_[2*idParticle + 1]\
                                         , xgl_[2*idParticle + 0], xgl_[2*idParticle + 1]);
    }
      fclose(fp);
//////////////////////////////////////////////
  sprintf(wallfile, "result/wall%05d.dat",iter/OUTGAP);
  if((hp = fopen(wallfile, "w")) == NULL){printf("FAILED TO OPEN 2nd FILE.\n"); exit(1);};
      fprintf(hp, "%15e %15e %15e %15e\n",0.000,0.000, L_[0],L_[1]);
      fclose(hp);
//////////////////////////////////////////////
    sprintf (filename2, "dataset/data%.4f.dat",RAMDA);
    if (iter==0){
      if((gp = fopen(filename2, "w")) == NULL){printf("FAILED TO OPEN 3rd FILE.\n"); exit(1);};
    }else{
      if((gp = fopen(filename2, "a")) == NULL){printf("FAILED TO OPEN 3rd FILE.\n"); exit(1);};
    }
    fprintf(gp,"%15e %15e %15e %15e %15e %15e %15e %15e %15e\n"\
              ,iter*DT,msd,(msd/(iter*DT*4.0)),SIGMA*SIGMA/2.0\
              ,(msd/(iter*DT*4.0))/D,L_[0],L_[1],Epx,alpha);
    fclose (gp);
//////////////////////////////////////////////
    

}

int main()
{
  Particle particle;

  for(int iter = 0; iter <= MAXITER; iter++){
    particle.Output(iter);

    particle.Update(iter);
  }
  return 0;
}
