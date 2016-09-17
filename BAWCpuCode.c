/*
	This is a main CPU code for Maxeler FPGA board. The goal here is to
	calculate the American style call/put option pricing based on Whaley model.

	For brevity, we will assume that all external factors, such as malloc, fopen, 
	perform smoothly.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "MaxSLiCInterface.h"
#include "Maxfiles.h"

// Assume risk-free interest rate is 0.03%
#define risk_free 0.0003 

// Data-array for Maxeler Engine
float *maxeler_call_price;
float *maxeler_put_price;
float *maxeler_rate;
float *maxeler_sigma_test;
float *maxeler_price_arr;
float *maxeler_strike_arr;
float *maxeler_ttm;


// Parameter for CPU Code
FILE* file_stream;
FILE *output_stream;


int *second;
int *delivery_date;

double *time_to_maturity;
double *underlying;
double *bid;
double *ask;
double *strike;
double *dynamic_sigma;
double *midprice;


enum option_type{
	call_option,
	put_option
};


/*
	This function will calculate the CDF of normal distribution ~ N(0,1).
	From prof. Bernt A Odegaard's book

	Parameter
		double *z: value to calculate for cdf

	Return
		double n: approximated cdf of normal distribution of z
*/
inline double N(double *z) {
    if ((*z) > 6.0) { return 1.0; }; // this guards against overflow
    if ((*z) < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs((*z));
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-(*z))*((*z)/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( (*z) < 0.0 ) n = 1.0 - n;
    return n;
}


/*
	This function will calculate the price of European call option
	Based on Black-Scholes-Merton model.

	Adapted from prof. Bernt A Odegaard's book

	Parameter
		const double *S:       // underlying spot price
		const double *K:       // strike price,
		const double *r:       // risk-free interest rate
		const double *sigma:   // volatility
		const double *time:    // time to maturity

	Return
		double n: Theoretical price of European call option
*/
inline double option_price_call_black_scholes (const double *S,      
                                               const double *K,    
                                               const double *r,      
                                               const double *sigma,  
                                               const double *time){ 
    double time_sqrt = sqrt(*time);
    double d1 = (log((*S)/(*K))+(*r)*(*time))/((*sigma)*time_sqrt)+0.5*(*sigma)*time_sqrt;
    double d2 = d1-((*sigma)*time_sqrt);
    double *p_d1 = &d1;
    double *p_d2 = &d2;
    return (*S)*N(p_d1) - (*K)*exp(-(*r)*(*time))*N(p_d2);
}


/*
	This function will calculate the price of European put option
	Based on Black-Scholes-Merton model.

	Adapted from prof. Bernt A Odegaard's book

	Parameter
		const double *S:       // underlying spot price
		const double *K:       // strike price,
		const double *r:       // risk-free interest rate
		const double *sigma:   // volatility
		const double *time:    // time to maturity

	Return
		double n: Theoretical price of European put option
*/
inline double option_price_put_black_scholes(const double *S,     
                                             const double *K,      
                                             const double *r,      
                                             const double *sigma,  
                                             const double *time){
    double time_sqrt = sqrt(*time);
    double d1 = (log((*S)/(*K))+(*r)*(*time))/((*sigma)*time_sqrt)+0.5*(*sigma)*time_sqrt;
    double d2 = d1-((*sigma)*time_sqrt);
    d1 = -d1;
    d2 = -d2;
    double *p_d1 = &d1;
    double *p_d2 = &d2;
    return (*K)*exp(-(*r)*(*time))*N(p_d2) - (*S)*N(p_d1);
}


/*

	This function will calculate the lagged underlying spot price based on the provided
	parameters. For example, with lag=-1 it will return the previous tick spot price.

	Parameter
		int current_t: Current tick number
		int lag: Desired lagged tick
		const double *underly: Pointer to underlying spot array
		const int *second: Pointer to trading data second array

	Return
		double p:	lagged spot price
*/
double get_lagged_underly(int current_t,int lag,const double *underly,const int *second){

	// if lag < 0, return previous tick data
    if(lag<0){
        if(current_t+lag>0)return underly[current_t+lag];
        else return -1;
    }

    // if lag = 0, return current tick data
    else if(lag == 0){
        return underly[current_t];
    }

    // if lag > 0, then keep looking for the lagged data
    else{

        int current_second = second[current_t];

        for(int i=current_t;i>0;i--){
            if(current_second-second[i]>=lag){
                return underly[i];
            }
        }

        // return -1 when lagged data not found
        return -1;

    }
}


/*

	This function will calculate the lagged strike price based on the provided
	parameters. For example, with lag=-1 it will return the previous tick strike price.

	Parameter
		int current_t: Current tick number
		int lag: Desired lagged tick
		const double *k: Pointer to strike price array
		const int *second: Pointer to trading data second array

	Return
		double p:	lagged strike price
*/
double get_lagged_strike(int current_t,int lag,const double *k,const int *second){

	
    if(lag<0){
        if(current_t+lag>0)return k[current_t+lag];
        else return -1;
    }

    else if(lag == 0){
        return k[current_t];
    }

    else{

        int current_second = second[current_t];

        for(int i=current_t;i>0;i--){
            if(current_second-second[i]>=lag){
                return k[i];
            }
        }
     
        return -1;

    }
}


/*

	This function will calculate the lagged time to maturity(TTM) based on the provided
	parameters. For example, with lag=-1 it will return the previous tick TTM.

	Parameter
		int current_t: Current tick number
		int lag: Desired lagged tick
		const double *t: Pointer to TTM array
		const int *second: Pointer to trading data second array

	Return
		double p:	lagged TTM
*/
double get_lagged_ttm(int current_t,int lag,const double *t,const int *second){

    if(lag<0){

        if(current_t+lag>0)return t[current_t+lag];
        else return -1;
    }

    else if(lag == 0){

        return t[current_t];
    }

    else{

        int current_second = second[current_t];

        for(int i=current_t;i>0;i--){
            if(current_second-second[i]>=lag){
                return t[i];
            }
        }

        return -1;

    }
}


/*
	This function will initilized the parameter for the main function, based
	on the type of option calculation.

	Parameter
		enum option_type type: type of option to calculate(call/put) 
		const char *filename_input: location of sample data
		const char *filename_output: location of desired output data
		const unsigned int size: size of sample data(number of lines)
*/
void initialized_parameter(enum option_type type, const char *filename_input,
 const char *filename_output, const unsigned int size){

	file_stream = fopen(filename_input, "r");

    output_stream = fopen(filename_output, "w");

    second = (int*)malloc(size*sizeof(int));

    time_to_maturity = (double*)malloc(size*sizeof(double));

    delivery_date=(int*)malloc(size*sizeof(int));

    underlying=(double*)malloc(size*sizeof(double));

    bid=(double*)malloc(size*sizeof(double));

    ask=(double*)malloc(size*sizeof(double));

    strike=(double*)malloc(size*sizeof(double));

    dynamic_sigma=(double*)malloc(size*sizeof(double));

    midprice=(double*)malloc(size*sizeof(double));

}


/*
	This function will import the necessary data from files and assign it
	into appropriate arrays for further operation.
*/
void import_data(){

	char *line = (char*)malloc(1024*sizeof(char));

    int i=0;
    //header = seconds,deliver,underlying,bid,ask,strike
    while (fgets(line, 1024, file_stream))
    {
        char* tok = strdup(line);
        for (tok = strtok(line, ",");
             tok && *tok;
             tok = strtok(NULL, ",\n"))
        {
            //second data
            sscanf(tok,"%d", &second[i]);
            tok = strtok(NULL, ",\n");

            //deliver_date
            int tmp_time;
            sscanf(tok,"%d", &tmp_time);
            delivery_date[i]=tmp_time;

            //TTM
            int year=(tmp_time-1503)/100;
            time_to_maturity[i]=(year*365+(tmp_time/100-year-15)*30+40)/365.0;
            tok = strtok(NULL, ",\n");

            //strike
            sscanf(tok,"%lf", &strike[i]);
            strike[i]/=100.0;
            tok = strtok(NULL, ",\n");

            //bid
            sscanf(tok,"%lf", &bid[i]);
            bid[i]/=100.0;
            tok = strtok(NULL, ",\n");

            //ask
            sscanf(tok,"%lf", &ask[i]);
            ask[i]/=100.0;
            tok = strtok(NULL, ",\n");

            //Underlying
            sscanf(tok,"%lf", &underlying[i]);
            underlying[i]/=100.0;
            tok = strtok(NULL, ",\n");

            //mid_price
            midprice[i]=(bid[i]+ask[i])/2.0;

        }

        i++;

        free(tok);
    }

    fclose(file_stream);

    printf("Read option data completed.\nStart reading Underlying asset data for implied volatility calculation.\n");

}


/*
	This function will perform the bisection method to calculate the implied volatility.

	Parameter
		unsigned int size: number of tick data
		enum option_type type: type of option to calculate(call/put)

*/
void calculate_implied_volatility(unsigned int size, enum option_type type){

	unsigned int idx = 0;

	double interest_rate = risk_free;

	// we will use bisection method to calculate volatility for each line
	while(idx<size)
	{
		// parameters setup
		int iteration_count = 0;

	    double low_vol = 0.001;
	    double high_vol = 1;
	    double epsilon = 0.000001;

	    // we assume that the target price is bid-ask midprice
	    double market_price = midprice[idx];
	    
	    // get necessary lagged data
	    double bst_underly = get_lagged_underly(idx,-1,underlying,second);
	    double bst_strike = get_lagged_strike(idx,0,strike,second);
	    double bst_ttm = get_lagged_ttm(idx,0,time_to_maturity,second);

	    // set underlying spot to current price, if it is the first tick
	    if(bst_underly==-1){
	        bst_underly=underlying[idx];
	    }

	    // set strike price to current price, if it is the first tick
	    if(bst_strike==-1){
	        bst_strike=strike[idx];
	    }

	    // set TTM to current value, if it is the first tick
	    if(bst_ttm==-1){
	        bst_ttm=time_to_maturity[idx];
	    }

	    // initial seed value for bisection
	    double seed_sigma =  0.5 * (low_vol + high_vol);

	    // initial option price to be used in bisection
	    double temp_price = 0;

	    if(type){
	    	temp_price = option_price_put_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm);
	    }

	    else
	    {
	    	temp_price = option_price_call_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm);
	    }

	     //= type==0?option_price_call_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm)\
	    :option_price_put_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm);//g(x);

	    // main loop for bisection method
	    while(fabs(temp_price - market_price)>epsilon) {

	    	// just in case the function is not converged
	        if(iteration_count>2500){
     			printf("bisection not converged, discard current line.\n");
     			break;
	        }

	        // adjust the lower bound up, if the price is too low
	        if (temp_price < market_price) {
	            low_vol = seed_sigma;
	        }

	        // adjust the upper bound down, if the price is too high
	        else if (temp_price > market_price) {
	            high_vol = seed_sigma;
	        }

	        // set new seed_sigma value
	        seed_sigma = 0.5 * (low_vol + high_vol);

	        // calculate new option price based on new volatility value(seed_sigma )
			if(type){
		    	temp_price = option_price_put_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm);
		    }

		    else
		    {
		    	temp_price = option_price_call_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm);
		    }

		    iteration_count++;
	        //temp_price = type==0?option_price_call_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm)\
	        :option_price_put_black_scholes(&bst_underly,&bst_strike,&interest_rate,&seed_sigma,&bst_ttm);
	        
	        //printf("error = %f\n",fabs(temp_price - market_price));

	    }

	    // store the volatility value in array
	    dynamic_sigma[idx] = seed_sigma;

	    idx++;

	}

	printf("Implied Volatility calculation completed.\n");

}


/*
	Helper function to prepare parameter array to stream into Maxeler's DFU engine/

	Parameter
		unsigned int size: number of tick data
		float interest_rate: risk-free interest rate
*/
void prepare_data_array_for_maxeler(unsigned int size, float interest_rate){

	maxeler_price_arr = (float*)malloc(size*sizeof(float));

    maxeler_strike_arr = (float*)malloc(size*sizeof(float));

    maxeler_ttm = (float*)malloc(size*sizeof(float));

    maxeler_rate = (float*)malloc(size*sizeof(float));

    maxeler_sigma_test = (float*)malloc(size*sizeof(float));

    maxeler_call_price = (float*)malloc(size*sizeof(float));

    maxeler_put_price = (float*)malloc(size*sizeof(float));

    for(unsigned int i=0; i<size; i++){

        maxeler_price_arr[i] = underlying[i];

        maxeler_strike_arr[i] = strike[i];

        maxeler_ttm[i]=time_to_maturity[i];

        maxeler_rate[i]=interest_rate;

        maxeler_sigma_test[i]=dynamic_sigma[i];
    }

    	memset(maxeler_call_price, 0, size);

    	memset(maxeler_put_price, 0, size);

}

/*
	Helper function to call the DFU engine to calculate the option pricing
	based on the parameter array that was prepared earlier.

	Parameter
		unsigned int size: number of tick data
		unsigned int sizeBytes: size in bytes of each data in array
*/
void call_maxeler(unsigned int size, unsigned int sizeBytes){

	BAW(size,maxeler_strike_arr,sizeBytes,maxeler_price_arr,sizeBytes,maxeler_rate,sizeBytes,
		maxeler_sigma_test,sizeBytes,maxeler_ttm,sizeBytes,maxeler_call_price,sizeBytes,
		maxeler_put_price,sizeBytes);

}


/*
	Helper function to store the calculated option price into CSV file format.

	Parameter
		enum option_type type: type of option(call/put)
		unsigned int size: number of result line
		float interest_rate: risk free interest rate
*/
void print_result(enum option_type type, unsigned int size, float interest_rate){

unsigned int correctly_price_sample = 0;

 if(type)
    {
        //header printing
        fprintf(output_stream,"S,K,r,sigma,TTM,Call,Bid,Ask,Percent_diff\n");
        for(int i=0;i<size;i++)
        {
        	// for overvalue price printing
            if(maxeler_call_price[i]>ask[i]){
                fprintf(output_stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n",underlying[i],strike[i],\
                	interest_rate,dynamic_sigma[i],time_to_maturity[i],maxeler_call_price[i],\
                	bid[i],ask[i],fabs(maxeler_call_price[i]-ask[i])/ask[i]*100);

            }

            // for undervalue price printing
            else if(maxeler_call_price[i]<bid[i]){
                fprintf(output_stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n", underlying[i],strike[i],\
                	interest_rate,dynamic_sigma[i],time_to_maturity[i],maxeler_call_price[i],\
                	bid[i],ask[i],fabs(maxeler_call_price[i]-bid[i])/bid[i]*100);
            }

            // for perfectly value price printing
            else{

                correctly_price_sample++;

                float mid=(bid[i]+ask[i])/2.0;

                fprintf(output_stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n", underlying[i],\
                	strike[i],interest_rate,dynamic_sigma[i],time_to_maturity[i],\
                	maxeler_call_price[i],bid[i],ask[i],fabs(maxeler_call_price[i]-mid)/mid*100);
            }
        }
    }

    else{
        //header printing
        fprintf(output_stream,"S,K,r,sigma,TTM,Put,Bid,Ask,Percent_diff\n");

        for(int i=0;i<size;i++)
        {
            if(maxeler_put_price[i]>ask[i])
            {
                fprintf(output_stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n", underlying[i],\
                	strike[i],interest_rate,dynamic_sigma[i],time_to_maturity[i],\
                	maxeler_put_price[i],bid[i],ask[i],fabs(maxeler_put_price[i]-ask[i])/ask[i]*100);
            }
            else if(maxeler_put_price[i]<bid[i])
            {
                fprintf(output_stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n", underlying[i],\
                	strike[i],interest_rate,dynamic_sigma[i],time_to_maturity[i],\
                	maxeler_put_price[i],bid[i],ask[i],fabs(maxeler_put_price[i]-bid[i])/bid[i]*100);
            }
            else
            {
                correctly_price_sample++;

                float mid=(bid[i]+ask[i])/2.0;

                fprintf(output_stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n", underlying[i],\
                	strike[i],interest_rate,dynamic_sigma[i],time_to_maturity[i],\
                	maxeler_put_price[i],bid[i],ask[i],fabs(maxeler_put_price[i]-mid)/mid*100);
            }
        }
    }

    fclose(output_stream);

    printf("Number of accurately price option = %d of %d.\n",correctly_price_sample, size);


}


/*

	Main function will control the overall flow of program. The main function
	will load the data samples from files and prepare them to be used in 
	implied volatility calculation and, ultimately, in Maxeler DFU engine.

*/
int main()
{
	enum option_type type = call_option;

	unsigned int file_line_count = 7081; // 9065 lines for our put sample-data

	const char *input_file = "/home/demo/Desktop/call_merge_output.csv";

	const char *output_file = "/home/demo/Desktop/call_merge_output_result.csv";

	long data_size_byte = file_line_count * sizeof(float);

    initialized_parameter(type, input_file, output_file, file_line_count);

    import_data();

    calculate_implied_volatility(file_line_count, type);

    prepare_data_array_for_maxeler(file_line_count, risk_free);

    call_maxeler(file_line_count, data_size_byte);

    print_result(type,file_line_count,risk_free);

    return 0;
}
