#include "mex.h"
//#include "stdio.h"

/*Function takes a vector x, a sparsity K, and cluster count C
and returns the appropriate support and cost*/

void bestKCsparse(double* xx, double* supp,int N, int K, int C)
{     
    //Compute sum of xx
    double sum_xx = 0.0;
    for (int ii=0;ii<N;ii++)
    {
        sum_xx = sum_xx + *(xx+ii);
    }
    
    int *****bkc_sol;
    bkc_sol=new int**** [N];
    for (int i=0; i<N; ++i)
    {
        bkc_sol[i]=new int*** [N];
        for (int j=0; j<N; ++j) 
        {
            bkc_sol[i][j]=new int** [K+1];
            for (int k=0; k<K+1; ++k)
            {
                bkc_sol[i][j][k] = new int* [C+1];
                for (int l=0;l<C+1;++l)
                {
                    bkc_sol[i][j][k][l] = new int [N];
                    for (int n=0;n<N;n++)
                    {
                        bkc_sol[i][j][k][l][n]=0;
                    }
                }
            }
        }
    }
    
    double ****bkc_cost;
    bkc_cost=new double*** [N];
    for (int i=0; i<N; ++i)
    {
        bkc_cost[i]=new double** [N];
        for (int j=0; j<N; ++j) 
        {
            bkc_cost[i][j]=new double* [K+1];
            for (int k=0; k<K+1; ++k)
            {
                bkc_cost[i][j][k] = new double [C+1];
                for (int l=0;l<C+1;++l)
                {
                    bkc_cost[i][j][k][l] = 10*sum_xx;
                }
            }
        }
    }
    
  int cc = 0;
  int jj=0;
  for (int kk=0;kk<=K;kk++)
  {
    for (int ii=0;ii<N;ii++)
    {
        for (int dd= 0;dd<=(N-1-ii);dd++)
        {
            jj=ii+dd;
            double sum = 0;
            for (int ss=ii;ss<=jj;ss++)
            {
                sum=sum+*(xx+ss);
            }
            bkc_cost[ii][jj][kk][cc]= sum;
        }
    }
  }
  

int kk=0;
for (int cc=0;cc<=C;cc++)
{
    for (int ii=0;ii<N;ii++)
    {
        for (int dd=0;dd<=(N-1-ii);dd++)
        {
            jj=ii+dd;
            double sum = 0;
            for (int ss=ii;ss<=jj;ss++)
            {
                sum=sum+*(xx+ss);
            }
            bkc_cost[ii][jj][kk][cc]= sum;
        }
    }
}




 //End Init

for (int cc=1;cc<=C;cc++)
{
    for (int kk=1;kk<=K;kk++)
    {
        for (int dd=0;dd<N;dd++)
        {
            for (int ii=0;ii<=N-dd-1;ii++)
            {
                jj=ii+dd;
                
                
                if (dd+1<=kk)
                {
                   
                    bkc_cost[ii][jj][kk][cc]=0;
                    for (int count=ii;count<=jj;count++)
                    {
                        
                        //printf("%d %d %d %d %d\n",ii,jj,kk,cc,count);
                        bkc_sol[ii][jj][kk][cc][count]=1;
                    }
                }
                
                else
                {
                    for (int cl=0;cl<=cc;cl++)
                    {
                        for (int kl=0;kl<=kk;kl++)
                        {
                            for (int ss=ii;ss<=jj-1;ss++)
                            {                                                                
                                double ccost = bkc_cost[ii][jj][kk][cc];
                                double lcost = bkc_cost[ii][ss][kl][cl];
                                double rcost = bkc_cost[ss+1][jj][kk-kl][cc-cl];
                                
                                if (ccost>lcost+rcost)
                                {
                                    bkc_cost[ii][jj][kk][cc] = lcost+rcost;
                                    //Loop through all count?
                                    for (int count=0;count<N;count++)
                                    {
                                        bkc_sol[ii][jj][kk][cc][count] = bkc_sol[ii][ss][kl][cl][count] + bkc_sol[ss+1][jj][kk-kl][cc-cl][count];
                                    }
                                }       
                            }
                        }
                      }
                    }
                }
            }
        }
    }

    for (int ii=0;ii<N;ii++)
    {
        *(supp+ii) = bkc_sol[0][N-1][K][C][ii];
    }

for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N; ++j) 
        {
            for (int k=0; k<K+1; ++k)
            {
                for (int l=0;l<C+1;++l)
                {
                   
                        delete [] bkc_sol[i][j][k][l];
                       
                }
                delete [] bkc_sol[i][j][k];
            }
            delete [] bkc_sol[i][j];
        }
        delete [] bkc_sol[i];
    }
    delete [] bkc_sol;
    bkc_sol = 0;

    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N; ++j) 
        {
            for (int k=0; k<K+1; ++k)
            {
                delete [] bkc_cost[i][j][k];
            }
            delete [] bkc_cost[i][j];
        }
        delete [] bkc_cost[i];
    }
    delete [] bkc_cost;
    bkc_cost = 0;    
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *xx, *supp;
  int N;
  int K, C;
  int mrows,ncols;

  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=3) 
    mexErrMsgTxt("Three inputs required: xx, K, C");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check to make sure the last two input arguments are scalar */
  if( !mxIsDouble(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
    mexErrMsgTxt("Input K must be a scalar.");
  }
  
  if( !mxIsDouble(prhs[2]) || mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
    mexErrMsgTxt("Input C must be a scalar.");
  }
  
  /*  get the scalar input K */
  K = mxGetScalar(prhs[1]);
  
  /*  get the scalar input C */
  C = mxGetScalar(prhs[2]);  
  
  /*  create a pointer to the input matrix y */
  xx = mxGetPr(prhs[0]);
  
  /*  get the dimensions of the matrix input xx */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if (mrows>ncols)
      N = mrows;
  else
      N=ncols;
  
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols,mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
   supp = mxGetPr(plhs[0]);

  /*  call the C subroutine */
   bestKCsparse(xx,supp,N,K,C);
  
}
