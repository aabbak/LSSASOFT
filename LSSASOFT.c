/*
 * This program is free software; you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software Foundation; 
 * either version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU General Public License for more details. You should have received a copy of 
 * the GNU General Public License along with this library; if not, write to the Free Software 
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * Authors    : R. Alpay ABBAK
 * Address    : Technical University of Konya, Geomatics Engineering Dept., Konya, TURKEY
 * E-Mail     : aabbak@selcuk.edu.tr
 *
 * Compilation of the program under Linux platform: 
 * g++ LSSASOFT.cpp matris.cpp -o LSSASOFT
 *
 * Running the program with our sample data: 
 * ./LSSASOFT test.dat -d3/0.1/20.1/28.1 -l
 * 
 * Modified: 23.11.2021, Konya           v1.0   
 */
#include<math.h>
#include"matris.h"
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<unistd.h>
#define MAX 1001
#define MX 90
#define NP 300
#define pi 3.1415926535898
#define OPTIONS "d:f:lwh"
matris SPEC(matris &P,matris &T,matris &F,matris &W,int NF,matris &DAT,int NDAT,int LT,matris &PER,int NPER);
double BASE(int i, double t, matris &DAT, int NDAT, int LT, matris &PER, int NPER);
void HELP();
int main(int argc, char *argv[])
{//***** DEFINATIONS ******************************************************************************//
	FILE *stream=NULL;			// input file pointer
	int 	 i=0;				// index of time series
	int 	 k=0;				// index of time series
	int	NF=0;				// number of data
	int     NK=0;				// number of known constituents
	int     LT=0;				// linear trend is exist or not (1=exist)
	int   NDAT=1;				// number of datum shifts
	int   NPER=0;				// number of forced periods
	int  wflag=0;				// weighted data flag
	char 	c;				// temprorary option
	double 	   N=0;				// length of data
	double 	  CL=0.0;			// critical percentage at 99% conf. level 
	double 	  Nq=0.0;			// Nyqustik Frequency
        double 	   p=0.0;			// the most dominant period
        double  smax=0.0;			// spectral value the most dominant period's 	
	matris     T(MAX);	 		// time vector
	matris     F(MAX);			// observation vector
	matris     W(MAX);	 		// weight vector
	matris     P(MAX);			// periods whose spectral values are desired
	matris     S(MAX);			// spectral values of the periods
	matris    DAT(MX);			// times of the datum biases 
	matris    PER(MX);			// forced periods
    	const char *ayrac ="/";			// discriminant symbol
//***** OPTIONS ANALYSIS *************************************************************************//
	while((c=getopt(argc,argv,OPTIONS))!=-1)
        switch(c)
        {
        	case 'd':
                	NDAT=atoi(strtok(optarg,ayrac));
			for(i=1;i<=NDAT;i++)
				DAT(i)=atof(strtok(NULL,ayrac));
                	break;
           	case 'f':
                	NPER=atoi(strtok(optarg,ayrac));
			for(i=1;i<=NPER;i++)
				PER(i)=atof(strtok(NULL,ayrac));
			break;
            	case 'l':
                	LT=1;
                	break;
             	case 'w':
                	wflag=1;
                	break;
             	case 'h':
                	HELP();
                	exit(EXIT_SUCCESS);
            	default:
                	HELP();
                	exit(EXIT_FAILURE);
        }
//***** READ FILE and DEFINE LIMITS **************************************************************//
	if((stream=fopen(argv[optind],"r")) == NULL)     	// it consists T,F(T)
    	{
		printf("%s data file can not be opened!\n",argv[optind]);
		exit(1);
    	}
	for(i=1;i<MAX;i++)
		W(i)=1.0;
    	i=1;
    	while(!feof(stream))
    	{
		if(wflag) fscanf(stream,"%lf%lf%lf\n",&T(i),&F(i),&W(i));
		else fscanf(stream,"%lf%lf\n",&T(i),&F(i));
		i++;
    	}
    	fclose(stream);
    	Nq=2.0*(T(2)-T(1));			// nyquistik frequency
	NF=i-1;					// the number of the data 						
	N=T(i-1);				// the length of data including the gaps
	DAT(1)=T(1);				// default datum bias  
	k=0;
//***** ESTIMATE PERIODS USING NEW ALGORITHM *****************************************************//
	do
	{
		smax=0.0;
		NK=NDAT+LT+2*NPER;		// total number of known constituents
		CL=100.0/(1.0/(pow(0.0001,(-2.0/(NF-NK-2)))-1)+1);
		for(i=1;i<=NP;i++)
			P(i)=(NP-1)/((NP-i)/(N/2)+(i-1)/Nq);
		S=SPEC(P,T,F,W,NF,DAT,NDAT,LT,PER,NPER);
		for(i=2;i<NP;i++)
		{
			if(S(i-1)<S(i) && S(i+1)<S(i) && smax<S(i))
			{
				smax=S(i); 
				p=P(i);
			}
		}
		if(smax<CL) exit(1); 	
		for(i=1;i<=NP;i++)
			P(i)=(NP-1)/((NP-i)/(p+1)+(i-1)/(p-1));
		S=SPEC(P,T,F,W,NF,DAT,NDAT,LT,PER,NPER);
		for(i=2;i<NP;i++)
		{
			if(S(i-1)<=S(i) && S(i+1)<=S(i) && smax<S(i))
			{
				smax=S(i); 
				p=P(i);
			}
		}
		if(smax>100.0) exit(1);	
		NPER++;
		k++;
		PER(NPER)=p;
		printf("%10.5lf %10.5lf\n",PER(NPER),smax);
	}while(smax>CL);
	return 0;
}
matris SPEC(matris &P,matris &T,matris &F,matris &W,int NF,matris &DAT,int NDAT,int LT,matris &PER,int NPER)
{//***** SPECTRAL ANALYSIS ************************************************************************//
	int 	     i=0;				// index of time series
	int 	     j=0;				// index of known constituents
	int 	     k=0;				// index of known constituents
	int         NK=0;				// number of known constituents
	double FNORM=0.0;				// quadratic norm of residual f (fTf)
	double  FUNC=0.0;				// Temrorary value 
	double OMEGA=0.0;				// radial frequency
        double  FCOS=0.0;				// f*cos(omega)
        double  FSIN=0.0;				// f*sin(omega)
        double    CC=0.0;				//
        double    CS=0.0;				//
        double    SS=0.0;				//
        double    WT=0.0;				// omega*t
        double COSWT=0.0;				// cosinus(omega*t)
        double SINWT=0.0;				// sinus(omega*t)
	double   UAU=0.0;				//
        double   UAV=0.0;				//
        double   VAV=0.0;				//
        double   DET=0.0;				// determinant of the matrix
	matris    S(MAX);				// spectral value of the periods
	matris     C(MX);				// amplitudes of known constituents
	matris     U(MX);				// 
	matris     V(MX);				// 
//***** CONSTRUCTION OF NORMAL EQUATIONS *********************************************************//
	NK=NDAT+LT+2*NPER;				// total number of known constituents
	matris   A(NK,NK);				// design matrix
	matris   Q(NK,NK);				// inversion of design matrix
	matris    B(NK+1);				// normal equation known vector
	for(i=1;i<=NF;i++)
	{
		for(j=1;j<=NK;j++)
		{
			FUNC=BASE(j,T(i),DAT,NDAT,LT,PER,NPER);
			B(j)=B(j)+FUNC*F(i)*W(i);
			for(k=j;k<=NK;k++)
				A(k-1,j-1)=A(k-1,j-1)+FUNC*BASE(k,T(i),DAT,NDAT,LT,PER,NPER);
		}
	}
	for(j=1;j<=NK;j++)
		for(k=j;k<=NK;k++)
			A(j-1,k-1)=A(k-1,j-1);
	Q=invch(A);
	A=Q;	
	for(i=1;i<=NK;i++)
		for(j=1;j<=NK;j++)
			C(i)=C(i)+A(i-1,j-1)*B(j);
	for(i=1;i<=NF;i++)
		for(j=1;j<=NK;j++)
			F(i)=F(i)-C(j)*BASE(j,T(i),DAT,NDAT,LT,PER,NPER);
	for(i=1;i<=NF;i++)
		FNORM=FNORM+F(i)*F(i)*W(i);
//***** COMPUTATION OF SPECTRUM ******************************************************************//
	for(i=1;i<=NP;i++)
	{
		OMEGA=2.0*pi/P(i);
         	FCOS=0.0;
         	FSIN=0.0;
         	CC=0.0;
         	CS=0.0;
         	SS=0.0;
		if(NK>0) 
		{
			for(j=1;j<=NK;j++)
			{
				U(j)=0.0;
				V(j)=0.0;
			}
		}
		for(j=1;j<=NF;j++)
		{
			WT=OMEGA*T(j);
			COSWT=cos(WT);
		        SINWT=sin(WT);
             		FCOS=FCOS+F(j)*COSWT;
             		FSIN=FSIN+F(j)*SINWT;
             		CC=CC+COSWT*COSWT;
             		CS=CS+COSWT*SINWT;
             		SS=SS+SINWT*SINWT;
			if(NK>0)
			{
				for(k=1;k<=NK;k++)
				{
					FUNC=BASE(k,T(j),DAT,NDAT,LT,PER,NPER);
                   			U(k)=U(k)+FUNC*COSWT;
                   			V(k)=V(k)+FUNC*SINWT;
				}
			}
		}
		UAU=0.0;
          	UAV=0.0;
          	VAV=0.0;
		if(NK>0)
		{
			for(j=1;j<=NK;j++)
			{
				for(k=1;k<=NK;k++)
				{
					UAU=UAU+U(j)*A(j-1,k-1)*U(k);
					UAV=UAV+U(j)*A(j-1,k-1)*V(k);
					VAV=VAV+V(j)*A(j-1,k-1)*V(k);
				}
			}
		}
		S(i)=0.0;
		DET=(CC-UAU)*(SS-VAV)-(CS-UAV)*(CS-UAV);
		S(i)=100.0*((SS-VAV)*FCOS*FCOS-2.0*(CS-UAV)*FCOS*FSIN+(CC-UAU)*FSIN*FSIN)/(DET*FNORM);
	}
	return S;
}
double BASE(int i, double t, matris &DAT, int NDAT, int LT, matris &PER, int NPER)
{//***** COMPUTATION OF FUNCTIONAL VALUES OF KNOWN CONSTITUENTS ***********************************//
	int ind=0;					//temprorary index
	double result=0.0;
	if(i<=NDAT)					// datum biases 
	{
		result=1.0;
		if((i==NDAT) && (t>=DAT(i)))
			return result;
		if((i<=NDAT) && (t>=DAT(i)) && (t<=DAT(i+1)))
			return result;
		result=0.0;
		return result;
	}
	if(i<=NDAT+LT) 					// linear trend
	{
		result=t;
		return result;
	}
	if(i<=NDAT+LT+2*NPER)				// forced period
	{
		ind=(i-NDAT-LT+1)/2;
		if(i-NDAT-LT!=ind*2)
		{
			result=cos(2.0*pi*t/PER(ind));
			return result;
		}
		result=sin(2.0*pi*t/PER(ind));
		return result;
	}
}
void HELP()
{//***** HELP *************************************************************************************//
    fprintf(stderr,"\nLSSASOFT estimates periodicities in time series by the LSSA method.\n\n");
    fprintf(stderr,"USAGE:\n");
    fprintf(stderr,"     LSSASOFT  data[file] -d[<value>] -f[<value>] -l -w -h\n\n");
    fprintf(stderr,"PARAMETERS:\n");
    fprintf(stderr,"     data:     includes time and variable. If it has the weighted data,\n");
    fprintf(stderr,"               third line must be the weight\n\n");
    fprintf(stderr,"OPTIONS:\n");
    fprintf(stderr,"     -d<value> number of datum biases, and times of the datum biases\n");
    fprintf(stderr,"               default: 1/0.0                               \n");
    fprintf(stderr,"     -f<value> number of forced periods, and their magnitudes\n\n");
    fprintf(stderr,"     -l        time series has a linear trend\n\n");
    fprintf(stderr,"     -w        time series has the weighted data\n\n");
    fprintf(stderr,"     -h        prints this help\n\n");
}
