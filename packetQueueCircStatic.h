//code from http://www.algoritmia.net/articles.php?key=colas#arraycir

#include <string.h>
//#define MAX_COLA 100000

typedef struct sPacket
{
    int id;
    int source;
    int destination;
    int deliveredOnTime;
    int whereAmI;
    int maxSteps;
    int numSteps;
    int hasPassedNode;
    int previous;
} packet;

typedef struct tcola
{
  int elems;
  int virtualElems;
  int entrada, salida;
  packet elementos[MAX_COLA];
} staticQueue;

void create(staticQueue * q)
{
  q->virtualElems = q->elems = q->salida = q->entrada = 0 ;
}

int empty(staticQueue * q)
{
  return (q->elems + q->virtualElems == 0);
}

int full(staticQueue *q)
{
  return (q->elems == MAX_COLA);
}

int numElems(staticQueue *q)
{
  return q->elems + q->virtualElems;
}

int numVirtualElems(staticQueue *q)
{
  return q->virtualElems;
}

//modif
void add(staticQueue *q, packet * elem)
{
    
    if(q->virtualElems == 0){
        if(full(q)){
            q->virtualElems++;
        }else{
            q->elems++;

            //q->elementos[q->entrada] = elem;
            memcpy(&q->elementos[q->entrada],elem,sizeof(packet));

            q->entrada ++;

            if(q->entrada == MAX_COLA){
                q->entrada = 0;
            }
        }
    }else{
        //if we have virtual elemnts we add new ones after the last
        q->virtualElems++;
    }
}

//modif
packet * get(staticQueue * q)
{
  void * returnElem = NULL;
  
  if(q->elems>0){
      q->elems--;      
      
      returnElem = &q->elementos[q->salida];
      
      q->salida++;
      if (q->salida == MAX_COLA)
        q->salida = 0;
  }else{
    if(q->virtualElems > 0){
        printf("Error: packet * get(staticQueue * q) -> trying to extract virtual elements\n");
    }
  }
  
  return returnElem;
}

packet * consultIndex(staticQueue * q, int index)
{
  void * returnElem = NULL;
  int i = 0;
  
  if(q->elems>index){

      returnElem = &q->elementos[(q->salida + index)%MAX_COLA];
  }else{
    if((q->elems + q->virtualElems) > index){
        printf("Error: packet * consultIndex(staticQueue * q, int index) -> trying to consult virtual elements\n");
    }
  }
  
  return returnElem;
}

//modif
packet * consult(staticQueue * q)
{
  void * returnElem = NULL;
  
  if(q->elems>0){
      returnElem = &q->elementos[q->salida];
  }else{
    if((q->virtualElems) > 0){
        printf("Error: packet * consult(staticQueue * q) -> trying to consult a virtual element\n");
    }  
  }
  
  return returnElem;
}

// 
// char outputString[1000];
// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
// {
//   struct tcola cola;
//   int elem;
//   
//     sprintf(outputString,"fprintf('FastGetNextRandomNode_Linear: all zeros in array\\n')");
//     mexEvalString(outputString);        
//   
//   crear(&cola);
//   if (vacia(&cola)) printf("\nCola vacia.");
//   if (llena(&cola)) printf("\nCola llena.");
//   encolar(&cola, 1);
//   desencolar(&cola, &elem);
//   return 0;
// }