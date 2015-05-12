#include "ldb.decl.h"
#include <cstdlib>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <algorithm>
#include <numeric>


typedef std::vector<int>::iterator IntIter;

int nChares;
int nRandMin;
int nRandMax;
CProxy_Main mainProxy;

class Main: public CBase_Main 
{
private:
  CProxy_IntArr arrProxy;
  double avg;
  int org_xor;
public:
  Main(CkArgMsg *m)
  { 
    CkPrintf("Argc: %d\n",m->argc);
    if (m->argc != 4)
    { 
      CkPrintf("Usage: ./charmrun ./mp1 +p4 chare-array-size min max\n");
      CkExit();
    }
    nChares  = atoi(m->argv[1]);
    nRandMin = atoi(m->argv[2]);
    nRandMax = atoi(m->argv[3]);
    mainProxy = thisProxy;
    arrProxy = CProxy_IntArr::ckNew(nChares);
    CkStartQD(CkCallback(CkReductionTarget(IntArr,PrintLB ), arrProxy));
  }

  void PrefixSumDone(int total)
  {
    CkPrintf("PrefixSum Max:%d\n",total);
    avg = double(total)/nChares;
    arrProxy.CreateBuffer(avg);
  }

  void BuffCreationDone()
  {
    int *pStart = new int[nChares];
    int *pEnd   = new int[nChares];
    for(int i=0; i<nChares; ++i)
    {
      pStart[i] = i*avg;
      pEnd[i]   = (i+1)*avg - 1;
    }
    arrProxy.LoadBalance(nChares,pStart,pEnd);
  }
 
  void xor_original(int xor_)
  {
    CkPrintf("Original XOR:%d\n",xor_);
    org_xor = xor_;
  }

  void done(int xor_)
  {
    CkPrintf("AfterLB XOR:%d\n",xor_);
    if( xor_ == org_xor)
     CkPrintf("Checksum match\n");
    else
     CkPrintf("Checksum doesn't match\n");
    CkExit();
  }
};

class IntArr: public CBase_IntArr
{
private:
  int nIntCount;
  int nPrefixSum; 
  int *pArray; //originalArray
  int nPrefixStage;
  int nPrefixStageCount;
  int *pValueBuf, *pFlagBuf;
  int nStartIndex,nEndIndex;
  int *pLBArray; //LoadBalancedArray
  double dAvg;
public:
  IntArr(CkMigrateMessage*) {}
  IntArr()
  {
    nPrefixStage = 0;
    nPrefixStageCount = ceil(log2 (nChares));
    //CkPrintf("StageCount:%d\n",nPrefixStageCount);
    pValueBuf = new int[nPrefixStageCount];
    pFlagBuf  = new int[nPrefixStageCount];
    memset(pFlagBuf, 0 ,nPrefixStageCount*sizeof(int));
    srand(time(NULL) + thisIndex);
    nIntCount = (rand()% (nRandMax - nRandMin)) + nRandMin;
  //push_kernel(push_input); 
    CkPrintf("Init count index[%d]:%d\n",thisIndex,nIntCount);
    pArray = new int[nIntCount];
    for(int i=0; i< nIntCount;++i)
    {
      pArray[i] = rand() % 100;
    }
    int xor_ = std::accumulate(pArray,pArray+nIntCount,0,std::bit_xor<int>());
    contribute(sizeof(int),&xor_,  CkReduction::bitvec_xor,
        CkCallback(CkReductionTarget(Main, xor_original), mainProxy));
 
    nPrefixSum = nIntCount;
    PrefixRun(nPrefixSum);
  }
  
  void PrefixRun(int value)
  {
    if(nPrefixStage >= nPrefixStageCount)
    {
      CkPrintf("PrefixSum[%d]:%d\n",thisIndex,nPrefixSum);
      delete[] pValueBuf;
      delete[] pFlagBuf;
      //delete[] pArray;
      contribute(sizeof(int),&nPrefixSum,  CkReduction::max_int,
        CkCallback(CkReductionTarget(Main, PrefixSumDone), mainProxy));
    }
    else {
      int sendIndex = thisIndex + (1<<nPrefixStage);
      if(sendIndex < nChares)
        thisProxy[sendIndex].PrefixRecVal(nPrefixStage,value); 
      
      if((thisIndex -(1<<nPrefixStage)) < 0) 
      {
        nPrefixStage++;
        PrefixRun(value);
      }
    }
  }
  
  void PrefixRecVal(int stage, int value)
  { 
    pFlagBuf[stage] = 1;
    pValueBuf[stage] = value;
    for(int i=nPrefixStage; i<nPrefixStageCount; ++i)
    {
      if(pFlagBuf[i] == 0) break;
      nPrefixSum+=pValueBuf[i];
      nPrefixStage++;
      PrefixRun(nPrefixSum);
    }
  }
  
  void LoadBalance(int chares, int *pStart, int *pEnd)
  {
   
    int nCurrStartIndex = nPrefixSum - nIntCount;
    int nCurrEndCount   = nPrefixSum -1;
    //go through the range an assign the owner
    //one more pass through and send to the correct owner
    std::vector<int> owningChare;
    for(int rank=nCurrStartIndex; rank<=nCurrEndCount;++rank)
    {  
       for(int i = 0;i<nChares;i++)
       {
        int min = pStart[i];
        int max = pEnd[i];
        if( rank >=min && rank<=max )
        {  
          owningChare.push_back(i);
          break;
        }
       } 
    }
    std::stringstream ss;
    for (IntIter itr = owningChare.begin(); itr!=owningChare.end();++itr)
       ss << *itr << ",";
    
    //std::cout << "LB Owner Assign:[" << thisIndex << "]" << ss.str() << std::endl;
    IntIter itr;
    int currChareIndex = owningChare[0];
    int currIndex = 0;
    int sameChareCount = 0;
    for(itr = owningChare.begin() + 1;itr != owningChare.end(); ++itr)
    {
      
      if(*itr == currChareIndex) 
      {
        sameChareCount++;
        if(itr+1 == owningChare.end())
        {
          /* this kind of senario: 0,0,0,1,1,*1
           */
          SendMsg(currChareIndex,sameChareCount+1,pArray+currIndex,nCurrStartIndex+currIndex);
          //thisProxy[currChareIndex].RecvArrMsg(sameChareCount+1,pArray+currIndex,nCurrStartIndex+currIndex);
        }
        else  continue;
      }
      else
      {
        /* 0,0,0,*1,1,1 */
        SendMsg(currChareIndex,sameChareCount+1,pArray+currIndex,nCurrStartIndex+currIndex); 
        //thisProxy[currChareIndex].RecvArrMsg(sameChareCount+1,pArray+currIndex,nCurrStartIndex+currIndex);
        currChareIndex = *itr;
        currIndex+= sameChareCount + 1;
        sameChareCount = 0;
         if(itr+1 == owningChare.end())
        { 
          /*0,0,0,1,1,1,*2 */
          SendMsg(currChareIndex,sameChareCount+1,pArray+currIndex,nCurrStartIndex+currIndex);
 
          //thisProxy[currChareIndex].RecvArrMsg(sameChareCount+1,pArray+currIndex,nCurrStartIndex+currIndex);
        }
      }  
    }
    std::string sNos("");
    char ch;
    std::stringstream s;
    for (int i = 0; i < nIntCount; i++) 
    { 
      s << pArray[i] << ",";
      //sNos += ch;
    }
    //std::cout << "Before LB[" << thisIndex << "]: " << s.str() << std::endl; 
    //CkPrintf("Before LB[%d] %s",thisIndex,sNos.c_str());
  }

  void PrintLB()
  { 
    std::string sNos("");
    char ch;
    std::stringstream s;
    for (int i = 0; i < nEndIndex - nStartIndex +1; i++) 
    { 
      s << pLBArray[i] << ",";
      
    }
    //CkPrintf("After LB[%d] %s",thisIndex,sNos.c_str());
    //std::cout << "After LB[" << thisIndex << "]: " << s.str() << std::endl;
    int xor_ = std::accumulate(pLBArray,pLBArray+nEndIndex - nStartIndex +1,0,std::bit_xor<int>());
    delete[] pLBArray;
    contribute(sizeof(int),&xor_,  CkReduction::bitvec_xor,
        CkCallback(CkReductionTarget(Main, done), mainProxy));
 }
  void SendMsg(int chareIndex,int size,int *data,int startIndex)
  {
    int *data_tmp = new int[size];
    memcpy(data_tmp,data,size*sizeof(int));
    thisProxy[chareIndex].RecvArrMsg(size,data_tmp,startIndex);
 
  }

  void CreateBuffer(double avg)
  {
    dAvg = avg;
    nStartIndex = thisIndex*avg;
    nEndIndex   = (thisIndex + 1)*avg - 1; //floor of both values
    CkPrintf("LB count[%d]: %d-%d\n",thisIndex,nStartIndex,nEndIndex);
    pLBArray = new int[nEndIndex - nStartIndex +1];
    contribute(CkCallback(CkReductionTarget(Main, BuffCreationDone), mainProxy));
  }

  void RecvArrMsg(int n, int *data, int startIndex)
  {
    int destIndex = startIndex - nStartIndex;
    memcpy(pLBArray+destIndex,data,n*sizeof(int));
    
  }
};
	

#include "ldb.def.h"

