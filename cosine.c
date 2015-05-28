#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <gmp.h>

#define EXPECTED_ARGC  5
#define OPTIONAL_ARGC  1
#define EVEN           0
#define ODD            1
#define COORDINATOR    0

typedef struct term {
    int number;                /*Thread identifier*/
    mpf_t result;              /*Result of the computation of this term*/  
    int computed;              /*Boolean used by the barrier to synchronize all the threads*/
    struct term *all_terms;    /*The array with all the thread arguments. Used by the barrier.*/
    int size;                  /*Size of the array. Used by the barrier.*/
    int print;                 /*To know if the programming is running with the 'd' parameter*/
    unsigned int turn;         /*To know the sequence of the threads to calculate power and factorial*/
    unsigned int n;            /*To calculate the term. Copies the value from the global*/
    int sgn;                   /*To calculate the term. Copies the value from the global*/
} Term;

/*Globals*/
char stop_condition;           /*Stop condition. Can be f or m*/  
mpf_t *precision;              /*Will be used to determine when to stop the program considering the value of stop_condition.*/
unsigned int counter;          /*Counts the number of turns. Everytime the threads leave the barrier, we have a turn*/
unsigned int n;                /*Used to know the terms. It's the n in the cosine series of the enunciation*/
int sgn;                       /*Oscilates between the values 1 and -1*/
mpf_t *x;                      /*precise value of x, in radians*/
mpf_t *cosx;                   /*The value of the cosine of x*/
mpz_t *factorial_array;        /*Array used to store the latest values. Have size 2, for n even, value is stored in factorial_array[0]. For odd n, in factorial_array[1]*/
mpf_t *power_array;            /*Same logic of factorial_array, but it is used to calculate de powers*/
unsigned int order;            /*Used to give an order of execution to the threads*/
unsigned int enter_power;      /*When allowed, this variable will allow a certain thread to calculate the power of its term*/
unsigned int enter_factorial;  /*When allowed, this variable will allow a certain thread to calculate the factorial of its term*/
unsigned int finish;           /*Knows when the program have to stop*/

/*Locks to guarantee mutual exclusion*/
pthread_mutex_t vlock;         /*Protects some global variables*/
pthread_mutex_t plock;         /*Protects power_array array*/
pthread_mutex_t flock;         /*Protects factorial_array array*/

/*Functions prototypes*/
/*************************************/
/*Input related*/
void args_qt(int);
void integer_argument(char*);
int do_first_argument(char**);
int do_second_argument(char**);
int do_third_argument(char**);
char *do_fifth_argument(int, char**);
/*Initialization/ending related*/
void initiate_terms(int, Term*, char*);
void free_terms(int, Term*);
void create_threads(int, pthread_t*, Term*);
void join_threads(int, pthread_t*);
void get_x(char*);
void get_precision(int);
void initiate_cosx();
void initiate_arrays();
void free_arrays();
void initiate_locks();
void destroy_locks();
/*cosine computation relation*/
void sequential_compute();
void *term_compute(void*);
void barrier(Term*, Term*);
void do_factorial2n(mpf_t, unsigned int);
mpz_t *mpz_fact2n(int);
void do_power2n(mpf_t, unsigned int);
mpf_t *mpf_pow2n(int);
void do_sign(mpf_t, int);
void do_absolute(mpf_t);
void assign_and_update(Term*);
void reset_computation_control_globals(unsigned int);
/*Print related*/
void print_thread_id(Term*);
void print_cosx(int);
/*************************************/

int main(int argc, char **argv)
{
    int q;                   /*It is the quantity of threads.*/
    int p;                   /*The integer passed through he command line to represent the precision*/
    char *xstr;              /*String containing the x value in radians*/
    char *behaviour;         /*Determine program's behaviour. Normal run(if eq NULL), printing additional information run(if eq 'd') or sequential run(if eq 's')*/
    pthread_t *term_threads; /*Contains all the 'q' threads that will be used to compute the cosine*/
    Term *thread_args;       /*The arguments for the threads: a struct term type for each thread*/

    /*Check the input and feed the program correctly*/
    args_qt(argc);
    /*Get command line arguments*/
    q = do_first_argument(argv);
    stop_condition = do_second_argument(argv);
    p = do_third_argument(argv);
    xstr = argv[4];
    behaviour = do_fifth_argument(argc, argv);

    /*Get precision and assigns to the global precision variable*/
    get_precision(p);

    /*Get x, in radians*/
    get_x(xstr);

    /*Initialize the cosine(x) with 0 value. Every new term computed is added to cosx*/
    initiate_cosx();

    /*Initialize globals*/
    counter = 0; n = 0; sgn = 1; order = 1;
    enter_power = 1; enter_factorial = 1; finish = 0;

    /*Initiate the arrays of size 2 that will be used to compute powers and factorials*/
    initiate_arrays();

    /*Initialize the locks*/
    initiate_locks();

    /*Normal run of 'd' run*/
    if((behaviour == NULL || *behaviour == 'd') && q > 1)
    {
        /*The thread 0 is a coordinator thread. The rest of the q threads are the term threads*/
        term_threads = malloc((q + 1) * sizeof(*term_threads));
        /*Thread arguments (term structs)*/
        thread_args = malloc((q + 1) * sizeof(*thread_args));

        /*Do computation*/    
        initiate_terms(q, thread_args, behaviour);
        create_threads(q, term_threads, thread_args);
        join_threads(q, term_threads);
    }
    /*'s' run*/
    else sequential_compute();

    /*Free all dynamic allocated memory*/
    if((behaviour == NULL || *behaviour == 'd') && q > 1)
    {
        /*Free thread and thread-related memory*/
        free_terms(q, thread_args);
        free(term_threads); term_threads = NULL;
        free(thread_args); thread_args = NULL;
    }
    /*Free memory of used by the power and factorial arrays*/
    free_arrays();
    /*Free memory used by the locks*/
    destroy_locks();
    /*Free memory of the remaining globals*/
    mpf_clear(*x); free(x); x = NULL;
    mpf_clear(*precision); free(precision); precision = NULL;
    mpf_clear(*cosx); free(cosx); cosx = NULL;
    return 0;
}

/*Calculate the cosine in the sequential mode*/
void sequential_compute()
{
    mpf_t result;                /*Will store the result of every loop iteration and will be reset after add this result to cos(x)*/

    /*Initialize result*/
    mpf_init(result);
    if(stop_condition == 'f')
    {
        mpf_t cosx_i;            /*Used to store the cos(x) value before the computations of the terms of this turn*/

        /*Initialize cosx_i*/
        mpf_init(cosx_i);

        for(n = 0; !finish; n++, sgn *= -1)
        {   
            /*Save the value of cos(x)*/
            mpf_set(cosx_i, *cosx);

            /*Set/reset the result variable*/
            mpf_set_ui(result, 1);

            /*Calculate the new value of cos(x)*/
            /*Do power*/
            do_power2n(result, n);

            /*Do Factorial*/
            do_factorial2n(result, n);

            /*Do sign*/
            do_sign(result, sgn);

            /*Calculate the absolute value of the difference of both cos(x)*/
            do_absolute(cosx_i);

            /*Check if the computation is finished*/
            if(mpf_cmp(cosx_i, *precision) < 0) finish = 1;

            /*Prints the cos(x)*/
            print_cosx(1);
        }
        /*Free cosx_i*/
        mpf_clear(cosx_i);
    }
    else
    {
        for(n = 0; !finish; n++, sgn *= -1)
        {   
            /*Set/reset the result variable*/
            mpf_set_ui(result, 1);

            /*Calculate the new value of cos(x)*/
            /*Do power*/
            do_power2n(result, n);

            /*Do Factorial*/
            do_factorial2n(result, n);

            /*Check if the computation of this term satisfies stop condition*/
            if(mpf_cmp(result, *precision) < 0) 
            {
                do_sign(result, sgn);
                /*Finished computing*/
                finish = 1;
            }

            /*Do sign*/
            do_sign(result, sgn);

            /*Prints the cos(x)*/
            print_cosx(1);
        }        
    }
    printf("\n\nComputation finished!\n\nTotal number of turns: %d\n\nFinal value of the cosine is:\n\n", n);
    print_cosx(1);

    /*Free result*/
    mpf_clear(result);
}

/*Threads to calculate a term of the cosine series*/
void *term_compute(void *args)
{
    Term *term = ((Term*) args);
   
    while(!finish) barrier(term->all_terms, term);

    if(term->number == COORDINATOR)
    {
        printf("\n\nComputation finished!\n\nTotal number of turns: %d\n\nFinal value of the cosine is:\n\n", counter);
        print_cosx(1);
    }

    return NULL;
}

void barrier(Term *all_terms, Term *term)
{
    /*Coordinator*/
    if(term->number == COORDINATOR)
    {
        mpf_t cosx_i;              /*Used to store the cos(x) value before the computations of the terms of this turn*/

        /*Print the threads as they reach this point of the program if running in mode d*/
        print_thread_id(term);

        /*Check if the stop condition is satisfied in case the program is running the 'f' mode*/
        if(stop_condition == 'f')
        {
            /*Initialize cosx_i*/
            mpf_init(cosx_i);
            /*Save the value of cos(x)*/
            mpf_set(cosx_i, *cosx);
        }

        term->all_terms[COORDINATOR].computed = 1;
        while(term->all_terms[term->size].computed != 1) { if(finish) return; continue; }

        /*Coordinator is counting the number of turns*/
        counter++;
        
        /*Check if the stop condition is satisfied in case the program is running the 'f' mode*/
        if(stop_condition == 'f')
        {
            /*Calculate the absolute value of the difference of both cos(x)*/
            do_absolute(cosx_i);
            /*Check if the computation is finished*/
            if(mpf_cmp(cosx_i, *precision) < 0) finish = 1;
            mpf_clear(cosx_i);
        }

        /*Prints the cosx*/
        print_cosx(term->print);

        /*'Reset' barrier to make it usable in the next iteration*/
        /*Used to guarantee synchronization of all the threads after the reset of the 'computed' field*/
        term->all_terms[COORDINATOR].sgn = 1;

        /*Print the threads as they reach this point of the program if running in mode d*/
        print_thread_id(term);

        /*Check if the stop condition is satisfied in case the program is running the 'f' mode*/
        if(stop_condition == 'f')
        {
            /*Initialize cosx_i*/
            mpf_init(cosx_i);
            /*Save the value of cos(x)*/
            mpf_set(cosx_i, *cosx);
        }

        term->all_terms[COORDINATOR].computed = 0;
        while(term->all_terms[term->size].computed != 0) { if(finish) return; continue; }

        /*Coordinator is counting the number of turns*/
        counter++;
        
        /*Check if the stop condition is satisfied in case the program is running the 'f' mode*/
        if(stop_condition == 'f')
        {
            /*Calculate the absolute value of the difference of both cos(x)*/
            do_absolute(cosx_i);
            /*Check if the computation is finished*/
            if(mpf_cmp(cosx_i, *precision) < 0) finish = 1;
            mpf_clear(cosx_i);
        }

        /*Prints the cosx*/
        print_cosx(term->print);

        /*Release all threads at the same time*/
        term->all_terms[COORDINATOR].sgn = 0;
    }

    /*Every other thread*/
    else
    {
        /*Print the threads as they reach this point of the program if running in mode d*/
        print_thread_id(term);
        
        /*As the threads grab this lock (inside the function), define the values that will be used to calculate this term and set the value for the next thread*/
        assign_and_update(term);

        /*Set/reset the result of this term*/
        mpf_set_ui(term->result, 1);

        /*Thread will wait its turn*/
        while(term->turn != enter_power) { if(finish) return; continue; }
        /*Do the power x^(2n) and add it (multiplicating) to the result*/
        do_power2n(term->result, term->n); enter_power++;
        
        /*Thread will wait its turn*/
        while(term->turn != enter_factorial) { if(finish) return; continue; }
        
        /*Do the factorial (2n)! and add it (dividing) to the result.*/
        /*Note: When the program reach this point, another thread will be calculating the power in parallel*/
        do_factorial2n(term->result, term->n); enter_factorial++;

        /*Wait the permission to execute*/
        while(term->all_terms[term->number - 1].computed == 0) { if(finish) return; continue; }

        /*Check if the stop condition is satisfied in case the program is running the 'm' mode*/
        if(stop_condition == 'm')
        {
            if(mpf_cmp(term->result, *precision) < 0) 
            {
                do_sign(term->result, term->sgn);
                /*Finished computing*/
                finish = 1;
            }
        }

        /*Add or sub the result to cosx, depending on the sgn*/
        do_sign(term->result, term->sgn);

        /*Last thread will reset globals*/
        reset_computation_control_globals(term->turn);

        /*Allow next term computation*/
        term->all_terms[term->number].computed = 1;
        
        /*Wait last term to be computed*/
        while(term->all_terms[term->size].computed != 1) { if(finish) return; continue; }

        /*Synchronize the thread id prints*/
        while(term->all_terms[COORDINATOR].sgn == 0) { if(finish) return; continue; }
        
        /*Print the threads as they reach this point of the program if running in mode d*/
        print_thread_id(term);
    
        /*As the threads grab this lock (inside the function), define the values that will be used to calculate this term and set the value for the next thread*/
        assign_and_update(term);

        /*Set/reset the result of this term*/
        mpf_set_ui(term->result, 1);

        /*Thread will wait its turn*/
        while(term->turn != enter_power) { if(finish) return; continue; }
        /*Do the power x^(2n) and add it (multiplicating) to the result*/
        do_power2n(term->result, term->n); enter_power++;
        
        /*Thread will wait its turn*/
        while(term->turn != enter_factorial) { if(finish) return; continue; }
        /*Do the factorial (2n)! and add it (dividing) to the result.*/
        /*Note: When the program reach this point, another thread will be calculating the power in parallel*/
        do_factorial2n(term->result, term->n); enter_factorial++;

        /*'Reset' barrier to make it usable in the next iteration*/
        while(term->all_terms[term->number - 1].computed == 1) { if(finish) return; continue; }

        /*Check if the stop condition is satisfied in case the program is running the 'm' mode*/
        if(stop_condition == 'm')
        {
            if(mpf_cmp(term->result, *precision) < 0) 
            {
                do_sign(term->result, term->sgn);
                /*Finished computing*/
                finish = 1;
            }
        }

        /*Add or sub the result to cosx, depending on the sgn*/
        do_sign(term->result, term->sgn);

        /*Last thread will reset globals*/
        reset_computation_control_globals(term->turn);

        term->all_terms[term->number].computed = 0;
        while(term->all_terms[term->size].computed != 0) { if(finish) return; continue; }
        while(term->all_terms[COORDINATOR].sgn == 1) { if(finish) return; continue; }
    }
}

/*Calculates the absolute value of the difference of the cos(x)*/
void do_absolute(mpf_t cosx_i)
{
    mpf_sub(cosx_i, *cosx, cosx_i);
    if(mpf_sgn(cosx_i) == -1) mpf_abs(cosx_i, cosx_i);
}

/*Resets the values of control global variables to the initial value of 1*/
void reset_computation_control_globals(unsigned int t)
{
    if(t == order)
    {
        order = 1; enter_factorial = 1; enter_power = 1;
    }
}

/*Prints the value of the cosx global variable*/
void print_cosx(int print)
{
    if(print == 1)
    {
        printf("\ncos(x) = ");
        mpf_out_str(stdout, 10, 0, *cosx);
        printf("\n-----------\n");
    }
}

/*Prints the thread id before it start a new computation*/
void print_thread_id(Term *term)
{
    if(term->print == 1 && !finish) printf("Thread%d ", term->number);
}

/*Do the sign and modify the cos(x) depending on it*/
void do_sign(mpf_t result, int s)
{
    if(!finish)
    {
        if(s == 1) mpf_add(*cosx, *cosx, result);
        else mpf_sub(*cosx, *cosx, result);
    }
}


/*Do the factorial of 2*i */
void do_factorial2n(mpf_t result, unsigned int i)
{   
    mpz_t *ptr_z;                /*Pointer used to retrieve the calculated factorial*/
    mpf_t temp;                  /*Used to convert mpz_t to mpf_t*/

    pthread_mutex_lock(&flock);
        mpf_init(temp);
            ptr_z = mpz_fact2n(2 * i);
            mpf_set_z(temp, *ptr_z);
            mpf_div(result, result, temp);
        mpf_clear(temp);
    pthread_mutex_unlock(&flock);
}

/*Assign values to the term use on calculation and prepare the global values to the next term to use*/
void assign_and_update(Term *term)
{
    pthread_mutex_lock(&vlock);
        /*Get the values from the globals*/
        term->sgn = sgn;
        term->n = n;
        term->turn = order;
        /*Set the variables for the next thread*/
        order++;
        n++;
        sgn *= -1;
    pthread_mutex_unlock(&vlock);
}

/*Do the power of 2*i*/
void do_power2n(mpf_t result, unsigned int i)
{   
    mpf_t *ptr_f;                /*Pointer used to retrieve the calculated power*/

    pthread_mutex_lock(&plock);
        ptr_f = mpf_pow2n(2 * i);
        mpf_mul(result, result, *ptr_f);
    pthread_mutex_unlock(&plock);
}


/*Factorial of 2n*/
mpz_t *mpz_fact2n(int tn)
{
    if(tn == 0) return &factorial_array[EVEN];
    else
    {
        mpz_mul_ui(factorial_array[EVEN], factorial_array[ODD], tn);
        mpz_mul_ui(factorial_array[ODD], factorial_array[EVEN], tn + 1);
        return &factorial_array[EVEN];
    }
}

/*Power of 2n*/
mpf_t *mpf_pow2n(int tn)
{
    if(tn == 0) return &power_array[EVEN];
    else
    {
        mpf_mul(power_array[ODD], power_array[EVEN], *x);
        mpf_mul(power_array[EVEN], power_array[ODD], *x);
        return &power_array[EVEN];
    }
}

/*Initiate arrays that will help calculate powers and factorials of 2n faster*/
void initiate_arrays()
{
    factorial_array = malloc(2 * sizeof(*factorial_array)); 
    mpz_init(factorial_array[EVEN]);
    mpz_set_ui(factorial_array[EVEN],1);
    mpz_init(factorial_array[ODD]);
    mpz_set_ui(factorial_array[ODD],1); 

    power_array = malloc(2 * sizeof(*power_array));
    mpf_init(power_array[EVEN]);
    mpf_set_ui(power_array[EVEN],1);
    mpf_init(power_array[ODD]);
    mpf_set(power_array[ODD], *x);
}

/*Free the memory allocated by the arrays used to help in the power and factorial calculations*/
void free_arrays()
{
    mpz_clear(factorial_array[EVEN]);
    mpz_clear(factorial_array[ODD]);
    free(factorial_array); factorial_array = NULL;
    
    mpf_clear(power_array[EVEN]);
    mpf_clear(power_array[ODD]);
    free(power_array); power_array = NULL;
}

/*Get the precision from the argument. Set it as a mpf_t type and set the default precision for the next calls of mpf_init*/
void get_precision(int p)
{
    unsigned long bits;
    char *string;
    char *mantissa = "1e";
    char exponent[1000], temp[1000];
    int flag;
    int c;

    /*Convert the number p (passed the commando line) to string*/
    string = malloc(1000 * sizeof(char));
    sprintf(temp, "%d", p);

    /*Build the string*/
    strcpy(exponent, "-");
    strcat(exponent, temp);
    strcat(strcpy(string, mantissa), exponent);


    /*For every 32bits, we gain ~10 more digits of precision*/
    /*c corrects the amount of digits we expect, because this relation is not precise*/
    c = p / 262;
    bits = ((32 * p)/ 10) + ((c + 1) * 32);

    /*Set the default precision for the next calls of mpf_init*/
    mpf_set_default_prec(bits);
    
    /*Initiate a mpf_t from the built string*/
    precision = malloc(sizeof(*precision));
    mpf_init(*precision);
    flag = mpf_set_str(*precision, string, 10);
    assert(flag == 0);

    free(string); string = NULL;
}

/*Initiate the cosx global*/
void initiate_cosx()
{
    cosx = malloc(sizeof(*cosx));
    mpf_init(*cosx);
    mpf_set_ui(*cosx, 0);
}

/*Convert and return the command line's 4th argument string to the x value, in radians*/
void get_x(char *xstr)
{
    float check;
    int flag;

    check = atof(xstr);
    if(check == 0)
    {
        x = malloc(sizeof(*x));
        mpf_init(*x);
        /*Seem precise enough to me.... This is 2pi*/
        flag = mpf_set_str(*x, "6.283185307179586476925286766559005768394338798750211641949889184615632812", 10);
        assert(flag == 0);   
    }
    else
    {
        x = malloc(sizeof(*x));
        mpf_init(*x);
        flag = mpf_set_str(*x, xstr, 10);
        assert(flag == 0);
    }
}

/*Free all the mpf_t type allocated in the threads*/
void free_terms(int q, Term *thread_args)
{
    int i;
    for(i = 0; i < q + 1; i++) mpf_clear(thread_args[i].result);
}

/*Initiate the terms, assigning a value to every field that needs a starting value or allocation*/
void initiate_terms(int q, Term *thread_args, char *behaviour)
{
    int i;

    thread_args[COORDINATOR].sgn = 0;
    for(i = 0; i < q + 1; i++) 
    { 
        thread_args[i].number = i;
        mpf_init(thread_args[i].result);
        thread_args[i].computed = 0;
        thread_args[i].size = q;
        thread_args[i].all_terms = thread_args;
        if(behaviour != NULL) thread_args[i].print = 1;
        else thread_args[i].print = 0;
    }
}

/*Function to create all term threads*/
void create_threads(int q, pthread_t *term_threads, Term *thread_args)
{
    int i;
    for(i = 0; i < q + 1; i++)
    {
        if (pthread_create(&term_threads[i], NULL, term_compute, &thread_args[i]))
        {
            printf("Error creating thread.");
            abort();
        }
    }
}

/*Function to join all term threads*/
void join_threads(int q, pthread_t *term_threads)
{
    int i;
    for(i = 0; i < q + 1; i++)
    {
        if (pthread_join(term_threads[i], NULL))
        {
            printf("Error joining thread.");
            abort();
        }
    }
}

/*Check if the number os arguments of the program are right*/
void args_qt(int argc)
{
    if(argc != EXPECTED_ARGC && (argc != EXPECTED_ARGC + OPTIONAL_ARGC))
    {
        printf("Wrong number of arguments. Use:\n\n\t./harmonic_cosine q [f|m] p x [d|s]\n\n\n");
        printf("Obs: last argument is optional.\n");
        exit(0);
    }
}

/*A fast check before calling atoi*/
void integer_argument(char *str)
{
    int i = 0;
    
    if(str[i] == '-') i++;
    if(str[i] == '\0')
    {
        printf("Error. Invalid input (first and third argument must be integers).\n");
        exit(0);   
    }
    
    while(str[i] != '\0')
    {
        if(str[i] >= '0' && str[i] <= '9') i++;
        else
        {
            printf("Error. Invalid input (first and third argument must be integers).\n");
            exit(0);       
        }
    }
}   

/*Initiate the 3 locks, used to guarantee mutual exclusion*/
void initiate_locks()
{
    if (pthread_mutex_init(&vlock, NULL) != 0)
    {
        printf("\nMUTEX initialization failed.\n");
        exit(1);
    }
    if (pthread_mutex_init(&plock, NULL) != 0)
    {
        printf("\nMUTEX initialization failed.\n");
        exit(1);
    }
    if (pthread_mutex_init(&flock, NULL) != 0)
    {
        printf("\nMUTEX initialization failed.\n");
        exit(1);
    }
}

/*Destroy the locks*/
void destroy_locks()
{
    pthread_mutex_destroy(&vlock);
    pthread_mutex_destroy(&plock);
    pthread_mutex_destroy(&flock);
}

int do_first_argument(char **argv)
{
    int q;
    integer_argument(argv[1]);
    if((q = atoi(argv[1])) < 0)
    {
        printf("Error. Found negative number passed as parameter (first argument).\n");
        exit(0);
    }
    if(q == 0) q = sysconf(_SC_NPROCESSORS_ONLN);
    return q;
}

int do_second_argument(char **argv)
{
    if(strlen(argv[2]) == 1 && (argv[2][0] == 'f' || argv[2][0] == 'm')) return argv[2][0];
    else
    {
        printf("Error. Found unexpected character. Was expecting f or m (second argument).\n");
        exit(0);
    }    
}

int do_third_argument(char **argv)
{
    int p;
    integer_argument(argv[3]);
    if((p = atoi(argv[3])) < 0)
    {
        printf("Error. Found negative number passed as parameter (third argument).\n");
        exit(0);
    }
    return p;
}

char *do_fifth_argument(int argc, char **argv)
{
    if(argc == EXPECTED_ARGC + OPTIONAL_ARGC)
    {
        if(strlen(argv[5]) == 1 && (argv[5][0] == 'd' || argv[5][0] == 's')) return &argv[5][0];
        else
        {
            printf("Error. Found unexpected character. Was expecting d or s (fifth argument).\n");
            exit(0);       
        }
    }
    return NULL;
}