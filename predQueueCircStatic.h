//code from http://www.algoritmia.net/articles.php?key=colas#arraycir

#include <string.h>
#include "mex.h"
#define TRUE 1
#define FALSE 0

typedef struct tcola
{
  int elems;
  int virtualElems;
  int entrada, salida;
  node_id_type elementos[MAX_COLA];
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
void add(staticQueue *q, node_id_type elem)
{
    
    if(q->virtualElems == 0){
        if(full(q)){
            q->virtualElems++;
        }else{
            q->elems++;

            //q->elementos[q->entrada] = elem;
            q->elementos[q->entrada] = elem;

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

node_id_type get(staticQueue * q)
{
  node_id_type returnElem = -1;
  
  if(q->elems>0){
      q->elems--;      
      
      returnElem = q->elementos[q->salida];
      
      q->salida++;
      if (q->salida == MAX_COLA)
        q->salida = 0;
  }else{
    if(q->virtualElems > 0){
        //printf("Error: packet * get(staticQueue * q) -> trying to extract virtual elements\n");
        mexErrMsgTxt("Error: packet * get(staticQueue * q) -> trying to extract virtual elements\n");
    }
  }
  
  return returnElem;
}


node_id_type pop(staticQueue * q)
{
  node_id_type returnElem = -1;
  
  if(q->elems>0){
      q->elems--;      
      
      q->entrada--;
      
      returnElem = q->elementos[q->entrada];      
      
      if (q->entrada == 0)
        q->entrada = MAX_COLA - 1;
  }else{
    if(q->virtualElems > 0){
        //printf("Error: packet * pop(staticQueue * q) -> trying to extract virtual elements\n");
        mexErrMsgTxt("Error: packet * pop(staticQueue * q) -> trying to extract virtual elements\n");
    }
  }
  
  return returnElem;
}


node_id_type consultIndex(staticQueue * q, int index)
{
  node_id_type returnElem = -1;
  int i = 0;
  
  if(q->elems>index){
      returnElem = q->elementos[(q->salida + index)%MAX_COLA];
  }else{
    if((q->elems + q->virtualElems) > index){
        //printf("Error: packet * consultIndex(staticQueue * q, int index) -> trying to consult virtual elements\n");
        mexErrMsgTxt("Error: packet * consultIndex(staticQueue * q, int index) -> trying to consult virtual elements\n");
    }
  }
  
  return returnElem;
}

//modif
node_id_type consult(staticQueue * q)
{
  node_id_type returnElem = -1;
  
  if(q->elems>0){
      returnElem = q->elementos[q->salida];
  }else{
    if((q->virtualElems) > 0){
        //printf("Error: packet * consult(staticQueue * q) -> trying to consult a virtual element\n");
        mexErrMsgTxt("Error: packet * consult(staticQueue * q) -> trying to consult a virtual element\n");
    }  
  }
  
  return returnElem;
}

int exist(staticQueue * q,node_id_type elem)
{
  int index = -1;
  int i = 0;
  if(q->elems>0){
      index = q->salida;
      while(i<q->elems){
        if(q->elementos[index] == elem ){
            return TRUE;
        }
        index ++;
        i ++;
      }
  
  }else{
    if((q->virtualElems) > 0){
        //printf("Error: packet * consult(staticQueue * q) -> trying to consult a virtual element\n");
        mexErrMsgTxt("Error: packet * consult(staticQueue * q) -> trying to consult a virtual element\n");
    }  
  }
  
  return FALSE;
}
 
// // to compile use mex -largeArrayDims testStaticQueue.c
// 
// #include "mex.h"
// #include "predQueueCircStatic.h"
// #include <string.h>
// #include <math.h>
// 
// char outputString[1000];
// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
//   struct tcola cola;
//   int elem;
//   
//     sprintf(outputString,"fprintf('FastGetNextRandomNode_Linear: all zeros in array\\n')");
//     mexEvalString(outputString);        
//   
//   create(&cola);
//   if (empty(&cola)) printf("\nCola vacia.\n");
//   if (full(&cola)) printf("\nCola llena.\n");
//   add(&cola, 1);
//   add(&cola, 2);
//   add(&cola, 3);
// 
//   sprintf(outputString,"fprintf('Get %i\\n')",get(&cola));
//   mexEvalString(outputString);        
//   sprintf(outputString,"fprintf('Pop %i\\n')",pop(&cola));
//   mexEvalString(outputString);        
//   sprintf(outputString,"fprintf('Pop %i\\n')",pop(&cola));
//   mexEvalString(outputString);        
//   sprintf(outputString,"fprintf('Get %i\\n')",get(&cola));
//   mexEvalString(outputString);          
//   sprintf(outputString,"fprintf('Pop %i\\n')",pop(&cola));
//   mexEvalString(outputString);        
//    
//   //desencolar(&cola, &elem);
// }