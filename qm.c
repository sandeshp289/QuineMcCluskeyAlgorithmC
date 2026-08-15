///////////////////////////////////////////////////////////////////////////////
// File       : qm.c
// Language   : C/C++
// Description: Quine-McCluskey algorithm
//		 - A step-by-step tabular method to simplify Boolean Functions
// Author     : Sandesh Prabhakar
// Created    : September 5, 2021
// Modified   : August 15, 2026
// Version    : 1.0
//
// Copyright (c) 2026 Sandesh Prabhakar
// All rights reserved
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct node
{
	int *value;
	node *link;
};

////////////////////////////////////////////////////////////////////////
// circuit class
////////////////////////////////////////////////////////////////////////

class circuit
{
    // circuit information
    int numvars;                 // total number of variables
    int minterms;                // total number of minterms
    int dontcares;               // total number of dontcares 
    int redun;                   // If 2, dont cares exists. If 1, dont cares do not exist
    int temp;                    // temporary variable to store minterms and dont cares while reading in
    int *minlist;    		 // minterm list
    int *redunlist;              // dont care list
    int **binlist;
    int **equalonecount;
    int *isdontcare;
    int ***table;
    int **tableM;
    int **index;
    int **compare;
    int *index_tableM;
    int combine_flag;
    int **MATRIX;
    int *MATRIX_COL;
    int **PI;
    int *MATRIX_ROW;
    int **EP;
    int *MATRIX_ROW_SORTED;
    int *MATRIX_ROW_COPY;
    int *epcount;

    // you may wish to add more private variables here

public:
    circuit(char *);   // constructor

    // you may wish to add more function declarations (or public variables) here
	int  decimal_to_binary (int** binlist, int* minlist, int minterms, int numvars, int** equalonecount, int*** table, int col_no, int** index );	
        int  compare_function (int** compare,int combine_flag,int numvars);
        int  sort_matrix_row(int* MATRIX_ROW_SORTED, int* MATRIX_ROW, int* MATRIX_ROW_COPY, int num_PI);
        void row_minimization(int* MATRIX_ROW_SORTED,int* MATRIX_ROW,int* MATRIX_ROW_COPY,int** MATRIX, int* MATRIX_COL, int num_PI,int minterms,int* minlist,int* epcount,int** PI, int** EP);

};

////////////////////////////////////////////////////////////////////////
// constructor: reads in the *.cubes file and builds basic data structure
// 		for the cube list 
////////////////////////////////////////////////////////////////////////
circuit::circuit(char *cktName)
{
    FILE *inFile;
    FILE *outFile;
    char fName[40];
    char foName[40];
    int i,j,k;
    char c;
    strcpy(fName, cktName);
    strcat(fName, ".sop");

    inFile = fopen(fName, "r");
    if (inFile == NULL)
    {
	fprintf(stderr, "Can't open .sop file\n");
	exit(-1);
    }

    fscanf(inFile, "%d", &numvars);	 // number of variables
    fscanf(inFile, "%d", &redun);    // if 2, this implies dontcares exist. If 1, there are no dont cares 

    int column_count,one_count,serial_number,table_content;
    column_count=numvars+10;
    one_count=numvars+1;
    serial_number=numvars*20;
    //printf("column count=%d one_count=%d serial_number=%d\n",column_count,one_count,serial_number);
    // allocate space for cubes data structure
	int max_cubes = int (pow (2,10)); 
        table_content=(serial_number+max_cubes)*5;

        minlist = new int [max_cubes+1];
        redunlist = new int [max_cubes+1];
        binlist = new int * [max_cubes+1];
        equalonecount = new int * [column_count+1];
        isdontcare = new int [max_cubes+1];
        table = new int ** [column_count+1];
        tableM = new int * [column_count+1];
        index = new int * [numvars+1];
        compare = new int * [numvars+1];
        index_tableM = new int [column_count+1];

        for(int i=0;i<=column_count;i++)
            {
            table[i]=new int * [one_count+1];
               for(int j=0; j<=one_count;j++)
                  {
                   table[i][j]=new int [table_content+1];
                  }
            }


          //Initializing the scheduling list
	 for(int i=0; i<=column_count;i++)
	 {
            for (int j=0; j<=one_count; j++)
	    {
               for (int k=0; k<=table_content; k++)
                    {
		    table[i][j][k] = 0;
                    }
            }
        }

	for (i=0; i<=max_cubes; i++)
	{
                binlist[i] = new int[numvars+2];
		for (j=0; j<=numvars+1; j++)
		{
		    binlist[i][j] = 0;
		}
		    redunlist[i] = 0;
		    minlist[i] = 0;
	}


	for (i=0; i<=column_count; i++)
	{
                tableM[i] = new int[serial_number+1];
		for (j=0; j<=serial_number; j++)
		    tableM[i][j] = 0;

                index_tableM[i]=0;
	}

	for (i=0; i<=numvars; i++)
	{
                compare[i] = new int[3];
		for (j=0; j<=2; j++)
		{
		    compare[i][j] = 0;
		}
	}

	for (i=0; i<=numvars; i++)
	{
                    index[i]=new int[column_count+1];
                    for(j=0;j<=column_count;j++)
		    index[i][j] = 0;
	}


	for (i=0; i<=column_count; i++)
	{
                equalonecount[i] = new int[numvars+1];
		for (j=0; j<=numvars; j++)
		{
		    equalonecount[i][j] = 0;
		}
	}


    // now read in the minterms
    fscanf(inFile, "%d", &minterms);     // number of minterms 
    for (i=1; i<=minterms; i++)
    {
	fscanf(inFile, "%d", &temp);
	minlist[i] = temp;
        isdontcare[temp]=0;
	}

    // now read in the dont cares if exist
	if(redun==2)
	{
  	  fscanf(inFile, "%d", &dontcares);     // number of minterms 

    		for (i=1; i<=dontcares; i++)
   		 {
		fscanf(inFile, "%d", &temp);
		redunlist[i] = temp;
                isdontcare[temp]=1;
		}

	}


	// skip till end of line
	while ((c = getc(inFile)) != '\n' && c != EOF)
	    ;
//    }	// for (i...)
    fclose(inFile);


    printf("Successfully read in minterms:\n");
    printf("\tPrinting minterm list\n");
    printf("\t%d No. of variables.\n", numvars);
    printf("\t%d No. of minterms\n", minterms);
    printf("\t");
    for (i=1; i<=minterms; i++)
    {
    printf("%d ", minlist[i]);
    }
    printf("\n");

    if(redun==2)
    {
    printf("\tPrinting dontcare list\n");
    printf("\t%d No. of dontcares\n", dontcares );
    printf("\t");
    for (i=1; i<=dontcares; i++)
    {
    printf("%d ", redunlist[i]);
    }
    printf("\n");
    }

    int col_no = 1;
    decimal_to_binary (binlist,minlist,minterms,numvars,equalonecount,table, col_no, index);


    printf("Printing binary values of the minterms:\n");
        for (i=1; i<=minterms; i++)
           {
           printf(" %d : ",minlist[i]);
            for (j=numvars-1; j>=0; j--)
             {
             printf("%d", binlist[minlist[i]][j]);
             }
             printf("\n");
           }


    if(redun==2)
    {
    decimal_to_binary (binlist,redunlist,dontcares,numvars, equalonecount,table, col_no, index);
    printf("Printing binary values of the dontcares:\n");
        for (i=1; i<=dontcares; i++)
           {
           printf(" %d : ",redunlist[i]);
            for (j=numvars-1; j>=0; j--)
             {
             printf("%d", binlist[redunlist[i]][j]);
             }
             printf("\n");
           }

    }


int flagT=0;
int mark=0;
int colT=col_no;
int z=0;
int duplicate_flag;

    while(!flagT)
     {
       mark=0;
       for(j=0;j<=numvars;j++)
             {
               z=0;
                while(z<index[col_no][j])
                  {
                  combine_flag=0;
                  int temp_flag=0;
                  z=z+colT; 
                    for(i=0;i<numvars;i++)
                      {
                      compare[0][i]=table[col_no][j][z+i];
                      }

                    if(j<numvars)
                     {

                     int k=0;
                     int qq;
                      while(k<index[col_no][j+1])
                        {
                  duplicate_flag=0;

                         k=k+colT;
                           for(qq=0;qq<numvars;qq++)
                             {
                             compare[1][qq]=table[col_no][j+1][k+qq];
                             }
                           combine_flag=compare_function(compare,combine_flag,numvars);

                           if(combine_flag==1)
                             {
                             temp_flag=1;

                             int x1=0;
                             int q1;
                             int tempq;
                             duplicate_flag=0;
                             tempq=pow(2,col_no);
                             while(x1<index[col_no+1][j])
                                 {
                                 duplicate_flag=0;
                                 x1=x1+tempq;
                                      for(q1=0;q1<numvars;q1++)
                                          {
                                          if(table[col_no+1][j][x1+q1]==compare[2][q1]) 
                                             {
                                             duplicate_flag++;
                                             }
                                          }
                                     if(duplicate_flag==numvars)
                                        {
                                        duplicate_flag=1;
                                        break;
                                        }
                                       else
                                        {
                                        duplicate_flag=0;
                                        }
                                 x1=x1+q1;
                                 }
                       if(duplicate_flag==0)
                          {
                              for(int r=0;r<colT;r++)
                                {
                                table[col_no+1][j][index[col_no+1][j]]=table[col_no][j][z-colT+r];
                                index[col_no+1][j]=index[col_no+1][j]+1;
                                }
                              for(int r=0;r<colT;r++)
                                {
                                table[col_no+1][j][index[col_no+1][j]]=table[col_no][j+1][k-colT+r];
                                index[col_no+1][j]=index[col_no+1][j]+1;
                                }
                              for(int r=0;r<numvars;r++)
                                {
                                table[col_no+1][j][index[col_no+1][j]]=compare[2][r];
                                index[col_no+1][j]=index[col_no+1][j]+1;
                                }
                         }
                             mark++;
                             }
                         k=k+qq;
                        }
                      }
                    if(j>0&&combine_flag==0&&temp_flag==0)
                       {
                        int k=0;
                        int qq;
                        while(k<index[col_no][j-1])
                         {
                         k=k+colT;
                           for(qq=0;qq<numvars;qq++)
                             {
                             compare[1][qq]=table[col_no][j-1][k+qq];
                             }
                           combine_flag=compare_function(compare,combine_flag,numvars);
                            if(combine_flag==1)
                            break;
                         k=k+qq;
                         }
                       }
                      if(combine_flag==0&&temp_flag==0)
                        {
                           int qq;
                           for(qq=0;qq<colT;qq++)
                           {
                           tableM[col_no][index_tableM[col_no]+qq]=table[col_no][j][z-colT+qq];
                           }
                           index_tableM[col_no]=index_tableM[col_no]+qq;

                           for(qq=0;qq<numvars;qq++)
                           {
                           tableM[col_no][index_tableM[col_no]+qq]=table[col_no][j][z+qq];
                           }
                           index_tableM[col_no]=index_tableM[col_no]+qq; 
                        }
                  z=z+i;
                  }

             }
     if(mark==0)
          {
          flagT=1;
          }
     else
         {
          col_no++;
          colT=colT*2;
         }
     }

printf("col_no = %d\n",col_no);
int col_count;
col_count=col_no;
col_no=1;
int y;
    printf("Listing all columns:\n");
 for(int n=1;n<=col_count;n++)
   {
    printf(" Column No. %d\n",n);
    printf(" No. of 1s        Decimal Label             Bit Strings\n");
       for (j=0; j<=numvars; j++)
           {
           int z=0; 
           while(z<index[n][j])
              {
              int tempx=pow(2,n-1);
              printf("   %d                   ",j);
              for(y=0;y<tempx;y++)
              {
              printf("%d,",table[n][j][z+y]);
              }
              z=z+y;
              printf("                  ");
                    for(i=numvars-1;i>=0;i--)
                      {
                        if(table[n][j][z+i]!=9)
                        {
                        printf("%d",table[n][j][z+i]);
                        }
                        else
                        {
                        printf("-");
                        }
                      }
                      printf("\n");
             z=z+numvars; 
              }
            printf("\n");
            }
   }


int num_PI=0;
    printf("Listing all prime implicants:\n");
    printf("Decimal Label             Bit Strings\n");
       for(int m=1;m<=col_count;m++)
            {
              int z=0;
              while(z<index_tableM[m])
                 {
                 int tempx=pow(2,m-1);
                     for(y=0;y<tempx;y++)
                      {
                      printf("%d,",tableM[m][z+y]);
                      }
                      z=z+y;
                      printf("                        ");
                    for(i=numvars-1;i>=0;i--)
                      {
                        if(tableM[m][z+i]!=9)
                        {
                        printf("%d",tableM[m][z+i]);
                        }
                        else
                        {
                        printf("-");
                        }
                      }
                      printf("\n");
                 z=z+numvars;
                 num_PI++;
                 }

             }

printf("No. of Prime Implicants=%d\n",num_PI);


MATRIX = new int * [max_cubes+1];
PI = new int * [num_PI+1];
EP = new int * [num_PI+1];
MATRIX_COL = new int [max_cubes+1];
MATRIX_ROW = new int [num_PI+1];
MATRIX_ROW_COPY = new int [num_PI+1];
MATRIX_ROW_SORTED = new int [num_PI+1];

epcount = new int [1];

epcount[0]=0;

	for (i=0; i<=max_cubes; i++)
	{
                MATRIX[i] = new int[num_PI+1];
		for (j=0; j<=num_PI; j++)
		{
		    MATRIX[i][j] = 0;
		}
		    MATRIX_COL[i] = 0;
	}

	for (i=0; i<=num_PI; i++)
	{
                PI[i] = new int[numvars+1];
                EP[i] = new int[numvars+1];
		for (j=0; j<=numvars; j++)
		{
		    PI[i][j] = 0;
		    EP[i][j] = 0;
		}
		    MATRIX_ROW[i] = 0;
		    MATRIX_ROW_COPY[i] = 0;
		    MATRIX_ROW_SORTED[i] = 0;
	}


int pi=0;
       for(int m=1;m<=col_count;m++)
            {
              int z=0;
              while(z<index_tableM[m])
                 {
                 pi++;
                 int tempx=pow(2,m-1);
                     for(y=0;y<tempx;y++)
                      {
                        int check=0;
                        for(int xy=1;xy<=dontcares;xy++)
                          {
                           if(tableM[m][z+y]==redunlist[xy])
                           check=1;
                          }
                          
                          if(check==0)
                          {
                          MATRIX[tableM[m][z+y]][pi]=1;
                          MATRIX_COL[tableM[m][z+y]]=MATRIX_COL[tableM[m][z+y]]+1;
                          MATRIX_ROW[pi]=MATRIX_ROW[pi]+1;
                          MATRIX_ROW_COPY[pi]=MATRIX_ROW_COPY[pi]+1;
                          }
                      }
                      z=z+y;
                      printf("                        ");
                    for(i=numvars-1;i>=0;i--)
                      {
                        if(tableM[m][z+i]!=9)
                        {
                        printf("%d",tableM[m][z+i]);
                        PI[pi][numvars-i]=tableM[m][z+i];
                        }
                        else
                        {
                        printf("-");
                        PI[pi][numvars-i]=tableM[m][z+i];
                        }
                      }
                      printf("\n");
                 z=z+numvars;
                 }

             }


             printf("\n          PRINTING INITIAL MATRIX\n");
             printf("          ");
       for(int m=1;m<=minterms;m++)
            {
             printf("%-6d    ",minlist[m]);
            }
             printf("\n");

       for(int m=1;m<=num_PI;m++)
            {
             for(i=1;i<=numvars;i++) 
              {
              if(PI[m][i]!=9)
              printf("%d",PI[m][i]);
              else
              printf("-");
              }
              printf("      ");
             for(j=1;j<=minterms;j++)
              {
              printf("%-6d    ",MATRIX[minlist[j]][m]);
              }
              printf("\n");
             }

      for(int m=1;m<=minterms;m++)
            {
                 if(MATRIX_COL[minlist[m]]==1)
                   {
                   for(i=1;i<=num_PI;i++)
                     {
                     if(MATRIX[minlist[m]][i]==1)
                       {
                       epcount[0]++;
                       for(j=1;j<=numvars;j++)
                          {
                          EP[epcount[0]][j]=PI[i][j];
                          }
                       MATRIX[minlist[m]][i]=0;
                       MATRIX_COL[minlist[m]]=0;
                       MATRIX_ROW[i]=MATRIX_ROW[i]-1;
                       MATRIX_ROW_COPY[i]=MATRIX_ROW_COPY[i]-1;
                         if(MATRIX_ROW[i]>0)
                           {
                            for(int yz=1;yz<=minterms;yz++)
                             {
                              if(MATRIX[minlist[yz]][i]==1)
                               {
                                MATRIX[minlist[yz]][i]=0;
                                MATRIX_COL[minlist[yz]]=MATRIX_COL[minlist[yz]]-1;
                                MATRIX_ROW[i]=MATRIX_ROW[i]-1;
                                MATRIX_ROW_COPY[i]=MATRIX_ROW_COPY[i]-1;

                                                        for(int p=1;p<=num_PI;p++)
                                                                 {
                                                                  if(MATRIX[minlist[yz]][p]==1)
                                                                         {
                                                                         MATRIX[minlist[yz]][p]=0;
                                                                         MATRIX_COL[minlist[yz]]=MATRIX_COL[minlist[yz]]-1;
                                                                         MATRIX_ROW[p]=MATRIX_ROW[p]-1;
                                                                         MATRIX_ROW_COPY[p]=MATRIX_ROW_COPY[p]-1;
                                                                         }
                                                                }
                               }
                             }
                           }
                       break;
                       } 
                     }
                   }
            }
             printf("\n");
             printf("\n          PRINTING REDUCED MATRIX\n");
             printf("          ");
       for(int m=1;m<=minterms;m++)
            {
             if(MATRIX_COL[minlist[m]]>0)
             printf("%-6d    ",minlist[m]);
            }
             printf("\n");

       for(int m=1;m<=num_PI;m++)
            {
            if(MATRIX_ROW[m]>0)
             {
             for(i=1;i<=numvars;i++) 
              {
              if(PI[m][i]!=9)
              printf("%d",PI[m][i]);
              else
              printf("-");
              }
              printf("      ");
             for(j=1;j<=minterms;j++)
              {
              if(MATRIX_COL[minlist[j]]>0)
              printf("%-6d    ",MATRIX[minlist[j]][m]);
              }
              printf("\n");
             }
            }


       for(i=1;i<=num_PI;i++)
            {
             MATRIX_ROW_SORTED[i]=i;
            }

            sort_matrix_row(MATRIX_ROW_SORTED,MATRIX_ROW, MATRIX_ROW_COPY,num_PI);

            row_minimization(MATRIX_ROW_SORTED,MATRIX_ROW,MATRIX_ROW_COPY,MATRIX, MATRIX_COL, num_PI,minterms,minlist,epcount,PI, EP);

             printf("Essential Primes\n");

       for(int m=1;m<=epcount[0];m++)
            {
             for(i=1;i<=numvars;i++)
              {
              if(EP[m][i]!=9)
              printf("%d",EP[m][i]);
              else
              printf("-");
              }
            printf("\n");
            }
             printf("\n");

}

////////////////////////////////////////////////////////////////////////
// other member functions here
////////////////////////////////////////////////////////////////////////
   int circuit::sort_matrix_row(int* MATRIX_ROW_SORTED, int* MATRIX_ROW, int* MATRIX_ROW_COPY, int num_PI)
        {
             int i,j,tempy, tempz;

             for(i=1;i<=num_PI;i++)
              {
               MATRIX_ROW_SORTED[i]=i;
              }

             for(i=1;i<=num_PI;i++)
              {
               MATRIX_ROW_COPY[i]=MATRIX_ROW[i];
              }

                for(i=1;i<num_PI;i++)
                 {
                   for(j=i+1;j<=num_PI;j++)
                     {
                       if(MATRIX_ROW_COPY[i]<MATRIX_ROW_COPY[j])
                         {
                         tempy= MATRIX_ROW_COPY[i];
                         MATRIX_ROW_COPY[i]=MATRIX_ROW_COPY[j];
                         tempz=MATRIX_ROW_SORTED[i];
                         MATRIX_ROW_SORTED[i]=MATRIX_ROW_SORTED[j];
                         MATRIX_ROW_COPY[j]=tempy;
                         MATRIX_ROW_SORTED[j]=tempz;
                        }
                    }
             }

        }



           void circuit::row_minimization(int* MATRIX_ROW_SORTED,int* MATRIX_ROW,int* MATRIX_ROW_COPY,int** MATRIX, int* MATRIX_COL, int num_PI,int minterms,int* minlist,int* epcount,int** PI, int** EP)
                {
                int i,j,k,temp; 
                temp=MATRIX_ROW_SORTED[1]; 
                epcount[0]++;
                for(i=1;i<=numvars;i++) 
                  {
                  EP[epcount[0]][i]=PI[temp][i]; 
                  }

                for(i=1;i<=minterms;i++)
                  {
                    if(MATRIX[minlist[i]][temp]==1)
                     {
                      MATRIX[minlist[i]][temp]=0;
                      MATRIX_COL[minlist[i]]=MATRIX_COL[minlist[i]]-1;
                      MATRIX_ROW[temp]=MATRIX_ROW[temp]-1;
                           for(j=1;j<=num_PI;j++)
                              {
                               if(MATRIX[minlist[i]][j]==1)
                                 {
                                 MATRIX[minlist[i]][j]=0;
                                 MATRIX_COL[minlist[i]]=MATRIX_COL[minlist[i]]-1;
                                 MATRIX_ROW[j]=MATRIX_ROW[j]-1; 
                                 }
                              }
                      }
                  }

                sort_matrix_row(MATRIX_ROW_SORTED,MATRIX_ROW, MATRIX_ROW_COPY,num_PI); 

             printf("\n");
               
             printf("\n");


                int empty=0;
                for(i=1;i<=num_PI;i++)
                 {
                  if(MATRIX_ROW_COPY[i]==0)
                   empty++;
                 }
                 if(empty!=num_PI)
                    {
                    row_minimization(MATRIX_ROW_SORTED,MATRIX_ROW,MATRIX_ROW_COPY,MATRIX, MATRIX_COL, num_PI,minterms,minlist,epcount,PI, EP);
                    }
                    else
                    {
                    return;
                    }

                }

		int circuit::decimal_to_binary (int** binlist, int* minlist, int minterms, int numvars, int** equalonecount, int*** table, int col_no, int** index )
                 {
				int i, j, z,temp, tempq;
                                int remainder;
                                int quotient;
                                int flag=1;
                                int onecount=0;
 
                           for (i=1; i<=minterms; i++)
                              {
                                 quotient=minlist[i];
                                 temp=quotient;
                                 if(temp!=0)
                                   {
                                    flag=0;
                                   }
                                   
                                for (j=0; j<numvars; j++) 
                                  {
                                    tempq=quotient;
                                    quotient=quotient/2;
                                    remainder=tempq%2;
                                    if (quotient != 0)
                                     {
                                     binlist[temp][j]=remainder; 
                                     if(remainder==1)
                                       {
                                       onecount++;
                                       }
                                     }
                                    else
                                     {
                                        if(flag==0)
                                          {
                                           binlist[temp][j]=1; 
                                           onecount++;
                                           flag=1;
                                          }
                                         else
                                          {
                                          binlist[temp][j]=0;
                                          }
                                     }
                                }
                                binlist[temp][numvars]=onecount;
                                equalonecount[col_no][onecount]++;
                                             table[col_no][onecount][index[col_no][onecount]]=minlist[i];
                                             index[col_no][onecount]++; 
                                           for (j=0; j<numvars;j++)
                                           {
                                               table[col_no][onecount][j+index[col_no][onecount]]=binlist[temp][j];
                                           }
                                           index[col_no][onecount]=j+index[col_no][onecount];
                                onecount=0;
                                            
                              }
                  }


     int circuit::compare_function(int** compare,int combine_flag,int numvars)
                     {
                     int count=0;
                     for(int i=0; i<numvars;i++)
                        {
                         if(compare[0][i]==compare[1][i])
                          {
                          compare[2][i]=compare[0][i];
                          }
                          else
                          {
                          compare[2][i]=9;
                          count++;
                          }
                        }
                     
                      if(count==1)
                          {
                          return 1;
                          }
                       else
                          {
                          return 0;
                          }
                     }

//
////////////////////////////////////////////////////////////////////////
// main starts here
////////////////////////////////////////////////////////////////////////
main(int argc, char *argv[])
{
    char cktName[20];
    circuit *ckt;
    int vecWidth, poNum;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <ckt>\n", argv[0]);
        exit(-1);
    }

    strcpy(cktName, argv[1]);

    // read in the cubes and build data structure
    ckt = new circuit(cktName);

    // NOW ADD YOUR FUNCTION HERE
}
