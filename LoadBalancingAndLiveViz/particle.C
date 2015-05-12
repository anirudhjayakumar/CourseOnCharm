#include <stdlib.h>
#include <vector>
#include "pup_stl.h"
#include "liveViz.h"
#include "Particle.h"
#include "ParticleExercise.decl.h"
#include <map>
#include <utility>
#define ITERATION (100)
#define MAX_INDEX  0
#define SUM_INDEX  1
#define ITER_INDEX 2
#define SIM_BOX_LEN 200

#define NW 0
#define N  1
#define NE 2
#define E  3
#define SE 4
#define S  5
#define SW 6
#define W  7



/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_Cell cellProxy;
/*readonly*/ int particlesPerCell;
/*readonly*/ int cellDimension;

//custom reduce typr
CkReduction::reducerType maxsum_IntType;
using namespace std;

class Main: public CBase_Main {
  private:
    double startTime;
  public:
    int doneCells;
    Main(CkArgMsg* m) {
      doneCells = 0;
      if(m->argc < 3) CkAbort("USAGE: ./charmrun +p<number_of_processors> ./particle <number of particles per cell> <size of array>");

      mainProxy = thisProxy;
      particlesPerCell = atoi(m->argv[1]);
      cellDimension = atoi(m->argv[2]);
      delete m;
      //CkPrintf("Creating grid of size %d\n",cellDimension); 
      CkArrayOptions opts(cellDimension,cellDimension);
      cellProxy = CProxy_Cell::ckNew(opts);
      
      //livevix init setup
      CkCallback cb_lv(CkIndex_Cell::getVizData(0),cellProxy);
      liveVizConfig cfg(liveVizConfig::pix_color,true);
      liveVizInit(cfg,cellProxy,cb_lv,opts);
      //CkPrintf("grid array created\n"); 
      cellProxy.run();
      startTime = CkWallTimer();
      
    }
    
    void done()
    {

      CkPrintf("Total time: %f\n",CkWallTimer() - startTime);
      CkExit();
    }

    //reductiontarget for the custom reduce
    void recvParticleCounts(CkReductionMsg *msg)
    {
      CkAssert(msg->getSize() == 3*sizeof(int));
      int *tmp = (int *) msg->getData();
      printTotal(tmp[SUM_INDEX],tmp[MAX_INDEX],tmp[ITER_INDEX]);
      delete msg;
    }
    // TODO: Add entry methods which will be a target of the reduction for avg
    // and max counts and exiting when the iterations are done

    void printTotal(int total, int max, int iter){
        CkPrintf("ITER %d, MAX: %d, TOTAL: %d\n", iter, max, total);
    }
};

// This class represent the cells of the simulation.
/// Each cell contains a vector of particle.
// On each time step, the cell perturbs the particles and moves them to neighboring cells as necessary.
class Cell: public CBase_Cell {
  Cell_SDAG_CODE
  private:
    int iteration;
    vector<Particle> particles;
    vector<Particle> NW_,N_,NE_,E_,SE_,S_,SW_,W_,L_; //L-lcoal
    map<int,pair<int,int> > mapNeighborToCord;
    map<pair<int,int>,int > mapCordToNeighbor;
    double minX,maxX,minY,maxY;
    int recv_count;
    unsigned char *pPixels;
    int pixelSize;
  public:
    Cell() {
      __sdag_init();
      iteration = 0;
      DefineBounds();
      CalcNeighbors();
      int squareSize = SIM_BOX_LEN/cellDimension;
      pixelSize = 3*squareSize*squareSize;
      pPixels = new unsigned char[pixelSize];
      usesAtSync = true;
      memset(pPixels,0,pixelSize);
      populateCell(particlesPerCell); //creates random particles within the cell
    }

    Cell(CkMigrateMessage* m) {}

    void pup(PUP::er &p){
      CBase_Cell::pup(p);
      __sdag_pup(p);
      p|iteration;
      p|particles;
      p| NW_;p|N_;p|NE_;p|E_;p|SE_;p|S_;p|SW_;p|W_;p|L_;
      p|mapNeighborToCord;
      p|mapCordToNeighbor;
      p|minX;p|maxX;p|minY;p|maxY;
      p|recv_count;
      p|pixelSize;
      if(p.isUnpacking())
      {
        pPixels = new unsigned char[pixelSize];
      }
      p(pPixels,pixelSize);
    }

    void UpdateParticles() {
      ClearTmpParticleHolders();
      for(PartIter iter = particles.begin();iter!=particles.end();++iter)
      {
        perturb(&(*iter));
        TransferParticleToTmpVector(&(*iter));
      }
      particles = L_;

    }

    void SendParticlesToNeighbors()
    {
      thisProxy(mapNeighborToCord[N].first,mapNeighborToCord[N].second).updateNeighbor(iteration,N_);
      thisProxy(mapNeighborToCord[E].first,mapNeighborToCord[E].second).updateNeighbor(iteration,E_);
      thisProxy(mapNeighborToCord[W].first,mapNeighborToCord[W].second).updateNeighbor(iteration,W_);
      thisProxy(mapNeighborToCord[S].first,mapNeighborToCord[S].second).updateNeighbor(iteration,S_);
      thisProxy(mapNeighborToCord[NW].first,mapNeighborToCord[NW].second).updateNeighbor(iteration,NW_);
      thisProxy(mapNeighborToCord[NE].first,mapNeighborToCord[NE].second).updateNeighbor(iteration,NE_);
      thisProxy(mapNeighborToCord[SW].first,mapNeighborToCord[SW].second).updateNeighbor(iteration,SW_);
      thisProxy(mapNeighborToCord[SE].first,mapNeighborToCord[SE].second).updateNeighbor(iteration,SE_);
    }
    
    
    void SendProgress()
    {
       int maxsum[3];
       maxsum[0] = maxsum[1] = particles.size();
       maxsum[2] = iteration;
       contribute(3*sizeof(int),maxsum,  maxsum_IntType, CkCallback(CkIndex_Main::recvParticleCounts(NULL), mainProxy));
    }

    void getVizData(liveVizRequestMsg *msg)
    {
      int squareSize = SIM_BOX_LEN/cellDimension;
      int locX = thisIndex.x*squareSize;
      int locY = thisIndex.y*squareSize;
      liveVizDeposit(msg,locX,locY,squareSize,squareSize,pPixels,this);
 
    }

      void ConstructFrame()
      {
        int squareSize = SIM_BOX_LEN/cellDimension;
               //CkPrintf("[%d,%d]SqSize:%d locX:%d locY:%d\n",thisIndex.x,thisIndex.y,squareSize,locX,locY);
        int totalPixels = squareSize*squareSize;
        
        memset(pPixels,0,3*totalPixels);
        for(PartIter iter = particles.begin();iter!=particles.end();++iter)
        {
          int partLocX = int(iter->x) % squareSize;
          int partLocY = int(iter->y) % squareSize;
          //CkPrintf("[%d,%d],(%d,%d)\n",thisIndex.x,thisIndex.y,partLocX,partLocY);
          if (iter->color == (int)e_BLUE)
          {
            pPixels[3*(partLocX*squareSize+partLocY)+0] = 0; // red
            pPixels[3*(partLocX*squareSize+partLocY)+1] = 0; // green
            pPixels[3*(partLocX*squareSize+partLocY)+2] = 255; // blue
          }
          else if (iter->color == (int)e_RED)
          {
            pPixels[3*(partLocX*squareSize+partLocY)+0] = 255; // red
            pPixels[3*(partLocX*squareSize+partLocY)+1] = 0; // green
            pPixels[3*(partLocX*squareSize+partLocY)+2] = 0; // blue
          }
          else 
          {
            pPixels[3*(partLocX*squareSize+partLocY)+0] = 0; // red
            pPixels[3*(partLocX*squareSize+partLocY)+1] = 255; // green
            pPixels[3*(partLocX*squareSize+partLocY)+2] = 0; // blue
          }
        }
      }
      
      /*
      // each cell will have diffent color pixels
      // proprtional to the particles 

      int totalBlue = GetTotalParticles(e_BLUE);
      int totalRed  = GetTotalParticles(e_RED);
      int totalGreen = GetTotalParticles(e_GREEN);

      int totalParticles = totalBlue + totalRed + totalGreen;
      CkPrintf("\nTotal particles: %d Green:%d Red:%d Blue:%d\n",totalParticles,totalGreen,totalRed,totalBlue);
      int bluePixels = (totalBlue/totalParticles) * totalPixels;
      int redPixels  = (totalRed/totalParticles) * totalPixels;
      int greenPixels = totalPixels - bluePixels - redPixels;


      usleep(100000);
      int counter = 0;
      for (int i=0;i<squareSize;++i)
      {
        for(int j=0; j<squareSize; ++j)
        {
          if (counter < bluePixels)
          {
            pPixels[3*(i*squareSize+j)+0] = 0; // red
            pPixels[3*(i*squareSize+j)+1] = 0; // green
            pPixels[3*(i*squareSize+j)+2] = 255; // blue
          }
          else if (counter < bluePixels + redPixels)
          {
            pPixels[3*(i*squareSize+j)+0] = 255; // red
            pPixels[3*(i*squareSize+j)+1] = 0; // green
            pPixels[3*(i*squareSize+j)+2] = 0; // blue
          }
          else 
          {
            pPixels[3*(i*squareSize+j)+0] = 0; // red
            pPixels[3*(i*squareSize+j)+1] = 255; // green
            pPixels[3*(i*squareSize+j)+2] = 0; // blue
          }
          counter++;

        }
      }
      */
            
  private:
        
    //private methods start here
    
    int GetTotalParticles(EColor eColor)
    {
      int counter = 0;
      for(PartIter iter = particles.begin();iter!=particles.end();++iter)
      {
        if (iter->color == (int)eColor) counter++;
      }
      return counter;
    }


    void populateCell(int initialElements) {
      //create random particles and add then to the particles vector
      if(thisIndex.x >= thisIndex.y) //blue zone
      {
        for (int i =0;i <initialElements;++i)
          particles.push_back(Particle(GenRandBtw(minX,maxX),GenRandBtw(minY,maxY),e_BLUE));
        //CkPrintf("[%d,%d] blue zone\n",thisIndex.x,thisIndex.y);
      }

      if(thisIndex.x <= thisIndex.y) //green zone
      {
        for (int i =0;i <initialElements;++i)
          particles.push_back(Particle(GenRandBtw(minX,maxX),GenRandBtw(minY,maxY),e_GREEN));
        //CkPrintf("[%d,%d] green zone\n",thisIndex.x,thisIndex.y);
      }

      //internel square dimension k
      int k = cellDimension/4;
      int startIndex = (cellDimension - k)/2;
      if( ((cellDimension-k)/2)%2 == 1) k++;
      int endIndex = startIndex + k -1;
      if (thisIndex.x >= startIndex && thisIndex.x <= endIndex \
          && thisIndex.y >= startIndex && thisIndex.y <= endIndex) { //red zone
        for (int i =0;i <2*initialElements;++i)
          particles.push_back(Particle(GenRandBtw(minX,maxX),GenRandBtw(minY,maxY),e_RED));
        //CkPrintf("[%d,%d] red zone\n",thisIndex.x,thisIndex.y);
      } 
      ConstructFrame();
      //CkPrintf("[%d,%d] Total particles:%d\n",thisIndex.x,thisIndex.y,particles.size());
    }



    /* wraping around the 2-d chare array */
    int wrapCell(int cell)
    {
      if(cell > (cellDimension -1 )) return cell - cellDimension;
      else if( cell < 0 ) return cell + cellDimension;
      else return cell;
    }

    /* wrapping around the box */
    double wrapBox(double point)
    {
      if(point > SIM_BOX_LEN)  return point - SIM_BOX_LEN;
      else if(point < 0.0)     return point + SIM_BOX_LEN;
      else return point;
    }

    /* transfers the particles to the transportation 
     * vectors called here as tmp vectors based on 
     * their cordinates 
     * */
    int TransferParticleToTmpVector(Particle *p)
    {
      double x = p->x;
      double y = p->y;
      x = wrapBox(x);
      y = wrapBox(y);
      EColor color = (EColor)p->color;
      double cellWidth = SIM_BOX_LEN/cellDimension;
      pair<int,int> cord(x/cellWidth,y/cellWidth);
      int dir = mapCordToNeighbor[cord];
      switch(dir)
      {
        case W:
          W_.push_back(Particle(x,y,color));
          break;
        case S:
          S_.push_back(Particle(x,y,color));
          break;
        case E:
          E_.push_back(Particle(x,y,color));
          break;
        case N:
          N_.push_back(Particle(x,y,color));
          break;
        case NW:
          NW_.push_back(Particle(x,y,color));
          break;
        case NE:
          NE_.push_back(Particle(x,y,color));
          break;
        case SW:
          SW_.push_back(Particle(x,y,color));
          break;
        case SE:
          SE_.push_back(Particle(x,y,color));
          break;
        default:
          L_.push_back(Particle(x,y,color));

      };
      //switch
      return 0;
    }

    void ClearTmpParticleHolders()
    {
      NW_.clear();N_.clear();NE_.clear();E_.clear();SE_.clear();S_.clear();SW_.clear();W_.clear();L_.clear(); 
    }

    // generate random alue between min and max
    double GenRandBtw(double min,double max)
    {
      return min + (drand48() * (max - min));
    }


    //change the location of the particle within the range of 8 neighbours
    //the location of the particles might exceed the bounds of the chare array
    //as a result of this functions, so you need to handle that case when deciding 
    //which particle to go which neighbour chare
    //e.g. the right neighbour of chare indexed[k-1,0] is chare [k,0]
    void perturb(Particle* particle) {
      //drand48 creates a random number between [0-1]	
      double deltax = (drand48()-drand48())*(SIM_BOX_LEN/cellDimension);
      double deltay = (drand48()-drand48())*(SIM_BOX_LEN/cellDimension);
      double speedFactor = 0;
  
      if(particle->color == e_RED) speedFactor = 1.0;
      else if (particle->color == e_GREEN) speedFactor = 0.25;
      else if (particle->color == e_BLUE) speedFactor = 0.50;
      else ReportErrorAndExit("Unknown color detected. Program terminating");

      particle->x += (deltax*speedFactor);
      particle->y += (deltay*speedFactor);
    }

    void ReportErrorAndExit(const char *msg)
    {
      CkPrintf(msg);
      CkExit();
    }
    //defines the boundry for each cell array element
    void DefineBounds()
    {
       minX = (SIM_BOX_LEN*thisIndex.x)/cellDimension;
       maxX = (SIM_BOX_LEN*(thisIndex.x+1))/cellDimension;
       minY = (SIM_BOX_LEN*thisIndex.y)/cellDimension;
       maxY = (SIM_BOX_LEN*(thisIndex.y+1))/cellDimension;
       //CkPrintf("DefineBounds[%d,%d]:X=[%f,%f],Y=[%f,%f]\n",thisIndex.x,thisIndex.y,minX,maxX,minY,maxY);
    }



    // each cell calculates its neighbors
    // this is a one-time operation during program 
    // startup
    void CalcNeighbors()
    {
      pair<int,int> coordinates;
      // west
      coordinates = std::make_pair(thisIndex.x,wrapCell(thisIndex.y-1));
      mapNeighborToCord[W] = coordinates;
      mapCordToNeighbor[coordinates] = W;
      //south
      coordinates = std::make_pair(wrapCell(thisIndex.x+1),thisIndex.y);
      mapNeighborToCord[S] = coordinates;
      mapCordToNeighbor[coordinates] = S;
      //east
      coordinates = std::make_pair(thisIndex.x,wrapCell(thisIndex.y+1));
      mapNeighborToCord[E] = coordinates;
      mapCordToNeighbor[coordinates] = E;
      //north
      coordinates = std::make_pair(wrapCell(thisIndex.x-1),thisIndex.y);
      mapNeighborToCord[N] = coordinates;
      mapCordToNeighbor[coordinates] = N;
      //nortwest
      coordinates = std::make_pair(wrapCell(thisIndex.x-1),wrapCell(thisIndex.y-1));
      mapNeighborToCord[NW] = coordinates;
      mapCordToNeighbor[coordinates] = NW;
      //northeast
      coordinates = std::make_pair(wrapCell(thisIndex.x-1),wrapCell(thisIndex.y+1));
      mapNeighborToCord[NE] = coordinates;
      mapCordToNeighbor[coordinates] = NE;
      //southwest
      coordinates = std::make_pair(wrapCell(thisIndex.x+1),wrapCell(thisIndex.y-1));
      mapNeighborToCord[SW] = coordinates;
      mapCordToNeighbor[coordinates] = SW;
      //southeast
      coordinates = std::make_pair(wrapCell(thisIndex.x+1),wrapCell(thisIndex.y+1));
      mapNeighborToCord[SE] = coordinates;
      mapCordToNeighbor[coordinates] = SE;

      /*
      typedef map<int,pair<int,int> >::iterator mapItr;
      for(mapItr it= mapNeighborToCord.begin(); it != mapNeighborToCord.end(); ++it)
      {
        CkPrintf("Neig->Coord[%d,%d]: %d -> [%d,%d]\n",thisIndex.x,thisIndex.y,it->first,it->second.first,it->second.second);
      }
      */
    }
};



/* custom reduction: each cell contributes three elements
 * index-0: max, index-1:sum,index-2:iteration
 * */
CkReductionMsg *maxsum_Int(int nmsg, CkReductionMsg **msgs)
{
  int maxsum[3]; 
  maxsum[MAX_INDEX] = 0;
  maxsum[SUM_INDEX] = 0;
  
  for(int j = 0; j < nmsg; j++){
    
    CkAssert(msgs[j]->getSize() == 3*sizeof(int));
    int *tmp = (int *) msgs[j]->getData();
    int value = tmp[MAX_INDEX];
    if(value > maxsum[MAX_INDEX]) maxsum[MAX_INDEX] = value;
    maxsum[SUM_INDEX]+=tmp[SUM_INDEX];
    if( j==0)
      maxsum[ITER_INDEX] = tmp[ITER_INDEX];
  }
 
 return CkReductionMsg::buildNew(3*sizeof(int), maxsum);
}
 
void register_maxsum_Int()
{
   maxsum_IntType = CkReduction::addReducer(maxsum_Int);
}



#include "ParticleExercise.def.h"
