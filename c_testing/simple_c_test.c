
#include <stdio.h>
#include <math.h>

int checkPrimeNumber(int n)
{
    int i;

    for(i=2; i <= n/2; ++i)
    {
        if(n%i == 0)
            return 1;
    }

    return 0;
}


int call_another_function(int n)
{
return n;
}