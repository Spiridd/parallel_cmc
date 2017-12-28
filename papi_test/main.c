#include <papi.h>
#include <stdlib.h>
#include <stdio.h>

#define NUM_EVENTS 2

void error_handler(const char* str)
{
    fprintf(stderr, "%s\n", str);
    exit(1);
}

int main()
{
    unsigned int Events[NUM_EVENTS] = 
            {PAPI_TOT_CYC, PAPI_BR_CN};
    int EventSet = PAPI_NULL;
    long long values[NUM_EVENTS];
    
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
        error_handler("cannot initialize");
    if (PAPI_create_eventset(&EventSet) != PAPI_OK)
        error_handler("cannot create eventset");
    if (PAPI_add_events(EventSet, Events, NUM_EVENTS) != PAPI_OK)
        error_handler("cannot add events");
    if (PAPI_start(EventSet) != PAPI_OK)
        error_handler("cannot start");

    // do work
    
    if (PAPI_stop(EventSet, values) != PAPI_OK)
        error_handler("cannot stop");
    
    printf("%lld PAPI_FP_INS\n", values[0]);
    printf("%lld PAPI_BR_CN\n", values[1]);

    const size_t EVENT_MAX = PAPI_num_counters();
    printf("#counters = %zu\n", EVENT_MAX);

    // is there branch counter?
    if (PAPI_query_event(PAPI_BR_CN) != PAPI_OK)
        printf("PAPI_BR_CN is not supported\n");
    // is there L2 data miss counter
    if (PAPI_query_event(PAPI_L2_DCM) != PAPI_OK)
        printf("L2 data cache misses are not supported\n");
    // is there L1 data miss counter
    if (PAPI_query_event(PAPI_L1_DCM) != PAPI_OK)
        printf("L1 data cache misses are not supported\n");
    
    return 0;
}

