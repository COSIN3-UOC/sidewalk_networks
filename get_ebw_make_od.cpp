//Code developed and used to produce the results of 
//
//  Rhoads, D., Solé-Ribalta, A., González, M. C., & Borge-Holthoefer, J. (2020). 
//  Planning for sustainable Open Streets in pandemic cities. 
//  arXiv preprint arXiv:2009.12548.
//
//Please cite the above paper in any works using or derived
//from this software


#define TRUE 1
#define FALSE 0
#define DEBUG 0
#include "mex.h"
#define MAX_COLA 65001
#define MAX_NODES 65000
typedef unsigned long int node_id_type;
#include "predQueueCircStatic.h"
#include <string.h>
#include <math.h>
#include "matrix.h"
#include <boost/heap/fibonacci_heap.hpp>

char outputString[150];

//returns the value of a half-normal PDF with sigma=855 at x
double halfNormalPDF(double x){

    //855 is the sigma that corresponds to a situation where
    //95% of trips fall within a distance of 1667 meters

    return (sqrt(2)/(855*(sqrt(M_PI)))) * exp(-pow(x,2)/(2*pow(855,2)));

}

//returns one row of the estimated OD matrix, 
//corresponding to the number of trips from s 
//to all other nodes j
void getOdRow(
    node_id_type s, 
    double * vStores, 
    double * vPop, 
    double * vCoords, 
    node_id_type sizeSupra, 
    double * attraction_vectPr, 
    double * OD_rowPr
){

    double xi, yi, xj, yj, dist, mass, rowSum;

    memset( attraction_vectPr, 0, sizeof(double)*sizeSupra );   
    memset( OD_rowPr, 0, sizeof(double)*sizeSupra );   

    //x and y coordinates of node s
    xi = vCoords[s];
    yi = vCoords[sizeSupra + s];

    //keep track of a sum of the row, to normalize later
    rowSum = 0;
    //for all nodes j in the network...
    for (node_id_type j=0; j<sizeSupra; j++){
        if (s != j ){

            //x and y coordinates of node j
            xj = vCoords[j];
            yj = vCoords[sizeSupra + j];

            //Calculate distance between s and j
            dist = sqrt(pow(xi-xj,2) + pow(yi-yj,2));

            //mass = # stores at j * halfnorm(distance_sj)
            mass = vStores[j] * halfNormalPDF(dist);
            attraction_vectPr[j] = mass;
            rowSum = rowSum + mass;

        }
        else {
            attraction_vectPr[j] = 0;
        }
    }

    for (node_id_type j=0; j<sizeSupra; j++){
        //OD_row = population at source * normalized mass of j
        OD_rowPr[j] = vPop[s]*(attraction_vectPr[j] / rowSum);
    }
}

void printMatrix(double * P, int rows, int cols);
void printMatrixInverted(double * P, int rows, int cols);
void printSparseMatrix(int n, double * pr, size_t * ir, size_t * jc);
//to use the heaps
struct dk_Node{
    node_id_type nodeId;
    double distance;
    dk_Node(node_id_type id, double dist) : nodeId(id),distance(dist) { } //<-- this is a constructor
};
struct compare_dk_Node{
    bool operator()(const dk_Node& n1, const dk_Node& n2) const{
        return n1.distance > n2.distance;
    }
};


//prhs[0] -> Adjacency Matrix in time
//prhs[1] -> N nodes
//prhs[2] -> N layers = 1
//prhs[3] -> vector s of services per node
//prhs[4] -> vector p of population per node
//prhs[5] -> vector c of node coordinates x,y

//plhs[0] -> Node betweenness
//plhs[1] -> Directed edge betweenness
//plhs[2] -> Num at start node CHECK
//plhs[3] -> Num at end node CHECK


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){


    //output matrices
    double tolerance = 0.000000001;
    node_id_type numNodes = 0, numLayers = 0;
    node_id_type sizeSupra = 0;
    mwSignedIndex dimsP[3];
    mxArray *rhs[1], *lhs[1];
    
    // get and check inputs
    if (nrhs != 6) mexErrMsgTxt( "Only 6 input arguments allowed." );
    if (nlhs > 2) mexErrMsgTxt( "Only 2 output argument allowed." );
    sizeSupra = mxGetN( prhs[0] );
    if (mxGetM( prhs[0] ) != sizeSupra) mexErrMsgTxt( "Input matrix G needs to be square." );
    if(mxIsSparse(prhs[0])==0) mexErrMsgTxt( "Input matrix must be sparse" );

    numNodes = mxGetPr(prhs[1])[0];
    numLayers = mxGetPr(prhs[2])[0];   

    double * vStores = mxGetPr(prhs[3]);
    double * vPop = mxGetPr(prhs[4]);
    double * vCoords = mxGetPr(prhs[5]);  

    mxArray * attraction_vect = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );
    mxArray * OD_row = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );
    double * attraction_vectPr = mxGetPr(attraction_vect);
    double * OD_rowPr = mxGetPr(OD_row);             
    
    if(numLayers != 1){
        mexErrMsgTxt("This code is not prepared to work in multiplex networks" );
    }
    if((numNodes*numLayers)!=sizeSupra) mexErrMsgTxt( "Incorrect size of the input matrix" );



    rhs[0] = mxDuplicateArray(prhs[0]);   
    mexCallMATLAB(1, lhs, 1, rhs, "'");
    mxDestroyArray(rhs[0]);

    
    double *Gpr = mxGetPr(lhs[0]);
    size_t *Gir = mxGetIr(lhs[0]);
    size_t *Gjc = mxGetJc(lhs[0]); 

    double edgeWeight = 0;
    double min_dist = -1;
    
    long int sIndEdges;
    long int eIndEdges;    
    long int nInd;
    node_id_type v,s,w,min_v;
    double numPacketReachDest;       
    double * C_B = 0;
    double * spEBW = 0;
    staticQueue * P = NULL;       
    double * sigma = 0;
    double * sigmaOrig = 0;    
    double * d = 0;
    double * delta = 0;
    double * visited = 0;
    staticQueue S;         
    
    //to store the nodes and the distance to them
    boost::heap::fibonacci_heap< dk_Node, boost::heap::compare<compare_dk_Node> > dk_heap;        
    boost::heap::fibonacci_heap< dk_Node, boost::heap::compare<compare_dk_Node> >::handle_type heapHandles[MAX_NODES];
    
    
    plhs[0] = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );
    C_B = mxGetPr(plhs[0]);    
    plhs[1] = mxCreateDoubleMatrix( numLayers*numNodes, numLayers*numNodes, mxREAL );
    spEBW = mxGetPr(plhs[1]);                      

    visited = (double *)mxMalloc(sizeof(double)*sizeSupra);     
    sigma = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    sigmaOrig = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    d = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    P = (staticQueue *)mxMalloc(sizeof(staticQueue)*sizeSupra);
    delta = (double *)mxMalloc(sizeof(double)*sizeSupra);     
    
    memset( C_B, 0, sizeof(double)*sizeSupra );   
    memset( spEBW, 0, sizeof(double)*sizeSupra );    
    
    for( s = 0 ; s<numNodes; s++){ 

        sprintf(outputString,"fprintf('Source = %lu\\n');",s+1);

        mexEvalString(outputString);

        getOdRow(s, vStores, vPop, vCoords, sizeSupra, attraction_vectPr, OD_rowPr);
        
        //initializations
        create(&S);  
              
        for( w = 0 ; w<sizeSupra ; w++){
            create(&P[w]);       
            sigma[w] = 0;                        
            d[w] = mxGetInf();
            visited[w] = 0;
    
            //store the handle in a table then use the table to access the handle
            heapHandles[w] = dk_heap.push(dk_Node(w,d[w]));             
        }            
        d[s] = 0;
        (*heapHandles[s]).distance = 0; dk_heap.decrease(heapHandles[s]);
        sigma[s] = 1;      
                 
        //end initializations  
                
        min_v = dk_heap.top().nodeId;
        min_dist = dk_heap.top().distance;
        dk_heap.pop();        
        
      
        

        while( min_dist < mxGetInf()){
        //Min_V := S;
        //Min_Dist := 0.0;
        //while Min_Dist /= Plus_Infinity loop

            if(DEBUG){
                sprintf(outputString,"fprintf('min_v %lu, min_dist = %f\\n');",min_v+1,min_dist);
                mexEvalString(outputString);                                                             
            }
            
            v = min_v;               
                                                                           
            
            visited[v] = 1;
            add(&S,v);                                                     

            
            if(DEBUG){
                sprintf(outputString,"fprintf('Exploring %lu\\n');",v+1);
                mexEvalString(outputString);                                                     
            }            
            
            //for each neigbor of v
            sIndEdges = Gjc[v];
            eIndEdges = Gjc[v+1] - 1;            

            if(DEBUG){
                sprintf(outputString,"fprintf('\\tmin_v = %lu, min_dist = %f\\n');",min_v+1, min_dist);
                mexEvalString(outputString);                                                                         
            }
                
            for( nInd = sIndEdges ; nInd<=eIndEdges ; nInd++ ){                        
            //while Has_Next(El) loop
                w = Gir[nInd];
                edgeWeight = Gpr[nInd];

                if(DEBUG){
                    sprintf(outputString,"fprintf('\tEdge from w=%lu to w=%lu and edgeWeight = %f\\n');",v+1, w+1, edgeWeight);
                    mexEvalString(outputString);                                                                         
                }                    
                
                if(DEBUG){                    
                    sprintf(outputString,"fprintf('\\tedges out %lu\\n');",w);
                    mexEvalString(outputString);     
                }

                if (v != w){
                //if V /= W then                    
                    if(d[w] > d[v] + edgeWeight){ // new path to the vertex
                        d[w] = d[v] + edgeWeight;
                        
                        (*heapHandles[w]).distance = d[w]; dk_heap.decrease(heapHandles[w]);                        

                        sigma[w] = sigma[v];
                        create(&P[w]);
                        add(&P[w],v);
                        
                        if(DEBUG){   
                            sprintf(outputString,"fprintf('\\t\\tnew path to vertex %lu with dist  = %f\\n');",w+1, d[w]);
                            mexEvalString(outputString);                                                     
                        }
                    }else if( fabs(d[w] - d[v] - edgeWeight) < tolerance){ // the current path is of equal length
                        sigma[w] = sigma[w] + sigma[v];
                        add(&P[w],v);                                                           
                    
                        if(DEBUG){                           
                            sprintf(outputString,"fprintf('\\t\\tpath to equal length to %lu , path difference  =%f\\n');",w+1, fabs(d[w] - d[v] - edgeWeight));
                            mexEvalString(outputString);          
                            sprintf(outputString,"fprintf('\\t\\told length = %f, dist to %lu = %f and edgeweight = %f\\n');",d[w],v+1,d[v],edgeWeight);
                            mexEvalString(outputString);                             
                        }
                    }                                                            
                 }
                 //end if;
            }
            //end loop;


            if(!dk_heap.empty()){
                min_v = dk_heap.top().nodeId;
                min_dist = dk_heap.top().distance;
                dk_heap.pop();                    
            }else{
                min_v = -1;
                min_dist = mxGetInf();            
            }
        }                        
                 
        //keep an original copy of the sigma
        //memcpy(sigmaOrig,sigma,sizeof(double)*sizeSupra);                
        for( w = 0 ; w<sizeSupra ; w++){
            delta[w] = 0;    
     
        }          
      
    
        while(!empty(&S)){
            
            w = pop(&S);
   

            if(DEBUG){
                sprintf(outputString,"fprintf('explore vertex w = %lu\\n');",w+1);
                mexEvalString(outputString);
            }
                        
            if(DEBUG){
                sprintf(outputString,"fprintf('Distance from (s=%lu) to (w=%lu) = %f\\n');",s+1,w+1,d[w]);
                mexEvalString(outputString);                                                         
            }
                
            while(!empty(&P[w])){ 

                v = pop(&P[w]);

            
                if(DEBUG){
                    sprintf(outputString,"fprintf('\\tvertex %lu has predecessor %lu\\n');",w+1,v+1);
                    mexEvalString(outputString);
                }
                
                if(DEBUG){
                    sprintf(outputString,"fprintf('\\t\\tPred: (v=%lu)\\n');",v);
                    mexEvalString(outputString);                                                         
                    sprintf(outputString,"fprintf('\\t\\tsigma[v=%lu]=%i/sigma[w=%lu]=%i,delta[w=%lu]=%f\\n');",v,(int)sigma[v],w,(int)sigma[w],w,delta[w]);
                    mexEvalString(outputString);                           
                }
                

                //--> Aquí la OD matrix


                numPacketReachDest = OD_rowPr[w];
  

                delta[v] = delta[v] + (sigma[v]/sigma[w])*(numPacketReachDest + delta[w]);  

                if(v != s){

                    spEBW[v + sizeSupra*w] = spEBW[v + sizeSupra*w] + (sigma[v]/sigma[w])*(numPacketReachDest + delta[w]);

                }                
                if(v == s){

                    spEBW[v + sizeSupra*w] = spEBW[v + sizeSupra*w] + (sigma[v]/sigma[w])*(numPacketReachDest + delta[w]);                                          
                  
                }                


                
            }


                
            if(w!=s){

                C_B[w] = C_B[w] + delta[w];  /// aqui la delta d'ha de normalizar no es poden sumar les betweenness
                                                               // de les differents capes simplement
                
                if(DEBUG){
                    sprintf(outputString,"fprintf('\\t\\tBetweenness[w(%lu)]=%f\\n');",w,C_B[w]);
                    mexEvalString(outputString);                                           
                }  
            
            }

        }

        
    }



    mxFree(sigma);
    mxFree(sigmaOrig);
    mxFree(d);
    mxFree(P);
    mxFree(delta);
    mxFree(visited);
}


void printMatrix(double * P, int rows, int cols){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < rows; rowIndex++){
        for(columnIndex = 0 ; columnIndex < cols; columnIndex++){            
            mexPrintf("%1.4f ",P[rowIndex*cols + columnIndex]);   
        }    
        mexPrintf("\n");   
    }        
}

void printMatrixInverted(double * P, int rows, int cols){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < rows; rowIndex++){
        for(columnIndex = 0 ; columnIndex < cols; columnIndex++){            
            mexPrintf("%1.4f ",P[columnIndex*rows + rowIndex]);   
        }    
        mexPrintf("\n");   
    }        
}

void printSparseMatrix(int n, double * pr, size_t * ir, size_t * jc){
    size_t nrow;

    mexPrintf("This is from the .cpp");

    for (int y=0; y<n; y++){
        nrow = jc[y+1] - jc[y];
        for (int x = 0; x<nrow; x++){
            mexPrintf("      (%d,%d)    %g\n", (*ir++)+1,y+1,*pr++);
        }
    }
}
    
