/**
 * Document: MaxCompiler Tutorial (maxcompiler-tutorial)
 * Chapter: 4      Exercise: 1      Name: Simple Exercise
 * MaxFile name: SimpleExercise
 * Summary:
 * 	 Takes a stream and for each value x calculates x^2 + x.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <MaxSLiCInterface.h>
#include <Maxfiles.h>


int main()
{
	FILE* stream = fopen("futuresTradeSorted.csv", "r");
	    int date[464130];
	    float time_to_maturity[464130];
	    float price[464130];
	    char line[1024];
	    int i=0;
	    while (fgets(line, 1024, stream))
	    {

	        char* tok = strdup(line);
	        for (tok = strtok(line, ",");
	             tok && *tok;
	             tok = strtok(NULL, ",\n"))
	        {
	            sscanf(tok,"%d", &date[i]);
	            tok = strtok(NULL, ",\n");
	            int tmp_time;
	            sscanf(tok,"%d", &tmp_time);
	            int year=(tmp_time-1503)/100;
	            time_to_maturity[i]=(year*365+(tmp_time/100-year-15)*30+40)/365.0;
	            tok = strtok(NULL, ",\n");
	            sscanf(tok,"%f", &price[i]);
	            tok = strtok(NULL, ",\n");
	        }
	        // NOTE strtok clobbers tmp
	        i++;
	        free(tok);
	    }

	    for (int i=0;i <100; i++)
	    {
	        printf("%d\t%f\t%f\t\n", date[i],time_to_maturity[i],price[i]);
	    }



    //price , strike, 0.01, 0, 0.28, 0.115
//	float para [6] = {100,105,0.01,0.0,0.28,0.115};
    float para[6]={0};
    para[0]=price[0]/100.0;
    para[1]=price[0]/100.0;
    para[2]=0.01;
    para[3]=0;
    para[4]=0.28;
    para[5]=time_to_maturity[0];
    printf("%f %f ",para[0],para[1]);
	float* y = malloc(sizeof(float)*8) ;

	Simple(para,64,y);

	printf("Call Price = %f\n", y[0]);

	return 0;
}