#include <iostream>
#include <cmath>
using namespace std;

const double G = 6.67384e-11;

struct body
{
double x, y, z;
double Vx, Vy, Vz;
double mass;
};

struct force
{
double Fx, Fy, Fz;
force():Fx(0), Fy(0), Fz(0){}
};

void input(body* data, int N);

int main(int argc, char* argv[])
{
   int dT, endtime, N;
   cin >> dT >> endtime >> N;
    
   body* data = new body[N];
   force* forces = new force[N];
   
   //DT is measured in seconds and is the length of one time tick.
   //endtime is the number of ticks to run the simulation
   //N is the number of bodies 
    
  
  
   //It's a start...need to finish declaring MPI_struct   
   MPI_Datatype MPI_Body;
   MPI_Datatype[7] =   {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    
   int blocklen[7] = {1, 1, 1, 1, 1, 1, 1};
   MPI_Aint disp[3]; 
    
   disp[0] = &partible
    
    
    
   input(data, N);
    
    //for each time tick
    for(int t = 0; t<endtime; t++)
    {
        //calculate the new forces
        for(int bOne = 0; bOne<N; bOne++)
        {
            force temp;
            for(int bTwo = 0; bTwo<N; bTwo++)
            {
                double dsquared = ((data[bOne].x-data[bTwo].x)*
                                   (data[bOne].x-data[bTwo].x)+
                                   (data[bOne].y-data[bTwo].y)*
                                   (data[bOne].y-data[bTwo].y)+
                                   (data[bOne].z-data[bTwo].z)*
                                   (data[bOne].z-data[bTwo].z));
                 
                if(dsquared != 0)
                {                  
                    double F = (G*data[bOne].mass*data[bTwo].mass)/dsquared;
                    
                    temp.Fx += F*(data[bTwo].x-data[bOne].x)/sqrt(dsquared);
                    temp.Fy += F*(data[bTwo].y-data[bOne].y)/sqrt(dsquared);
                    temp.Fz += F*(data[bTwo].z-data[bOne].z)/sqrt(dsquared);
                }
            }
            forces[bOne] = temp;
        }
        
        //calculate new velocities and locations
        for(int currentBody = 0; currentBody < N; currentBody++)
        {
            data[currentBody].Vx += forces[currentBody].Fx * dT / 
                                                        data[currentBody].mass;
            data[currentBody].Vy += forces[currentBody].Fy * dT / 
                                                        data[currentBody].mass;
            data[currentBody].Vz += forces[currentBody].Fz * dT / 
                                                        data[currentBody].mass;
            
            data[currentBody].x += data[currentBody].Vx * dT;
            data[currentBody].y += data[currentBody].Vy * dT;
            data[currentBody].z += data[currentBody].Vz * dT;
        }
        
        //output current locations!
        //commented out non-location output since this allows gnuplot plotting.
        if(t%30 == 0)
        {
            for(int currentBody = 0; currentBody < N; currentBody++)
            {
                cout //<< "body " << currentBody << ": " 
                     << data[currentBody].x << " "
                     << data[currentBody].y << " "
                     << data[currentBody].z << "  " << endl;
                    /* << data[currentBody].Vx << " "
                     << data[currentBody].Vy << " "
                     << data[currentBody].Vz << endl; */
            }        
        }
    }  
    
    delete[] data;
    delete[] forces;
    return 0;
}

void input(body* data, int N)
{
    body temp;
    for(int x = 0; x<N; x++)
    {
    cin >> temp.x >> temp.y >> temp.z >> 
          temp.Vx >> temp.Vy >> temp.Vz >> temp.mass;
    
        data[x] = temp;
    }
}


