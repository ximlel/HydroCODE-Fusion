/**
 * @file  terminal_io.c
 * @brief This is a set of common functions for terminal input and output.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/var_struc.h"


/**
 * @brief This is a functions preprocesses ARGuments.
 * @details This function prints out all ARGuments, checks for the right ARGument Counter, loads argv[3], argv[4] and
 *          argv[argc_least+1,argc_least+2,…] as configuration, and puts the scheme name into the pointer 'scheme'.
 * @param[in] argc_least: The least value of the ARGument Counter.
 * @param[in] argc:       ARGument Counter.
 * @param[in] argv:       ARGument Values.
 * @param[in] scheme:     Scheme name.
 *          - argv[1]: Folder name of test example (input path).
 *          - argv[3]: Dimensionality (= 1).
 *          - argv[4]: Order of numerical scheme[_scheme name] (= 1[_Riemann_exact] or 2[_GRP]).
 *          - argv[argc_least+1,argc_least+2,…]: Configuration supplement config[n]=(double)C (= n=C).
 */
void arg_preprocess(const int argc_least, const int argc, char *argv[], char * scheme)
{
    int k, j;
    printf("\n");
    for (k = 0; k < argc; k++)
	printf("%s ", argv[k]);
    printf("\n\n");
#ifdef _WIN32
    printf("TEST:\n %s\n", argv[1]);
#elif __linux__
    printf("\x1b[47;34mTEST:\x1b[0m\n \x1b[1;31m%s\x1b[0m\n", argv[1]);
#endif
    if(argc < argc_least)
	{
#ifdef _WIN32
	    printf("Test Beginning: ARGuments Counter %d is less than %d\n", argc, argc_least);
#elif __linux__
	    printf("Test Beginning: \x1b[43;37mARGuments Counter %d is less than %d\x1b[0m\n", argc, argc_least);
#endif
	    exit(4);
	}
    else
#ifdef _WIN32
	printf("Test Beginning: ARGuments Counter = %d\n", argc);
#elif __linux__
    printf("Test Beginning: \x1b[43;37mARGuments Counter = %d\x1b[0m\n", argc);
#endif

    // Set dimension.
    int dim;
    dim = atoi(argv[3]);
    config[0] = (double)dim;

    // Set order and scheme.
    int order; // 1, 2
#ifdef _WIN32
    printf("Order[_Scheme]: %s\n",argv[4]);
#elif __linux__
    printf("Order[_Scheme]: \x1b[41;37m%s\x1b[0m\n",argv[4]);
#endif
    errno = 0;
    order = strtoul(argv[4], &scheme, 10);
    if (*scheme == '_')
	scheme++;
    else if (*scheme != '\0' || errno == ERANGE)
	{
	    printf("No order or Wrog scheme!\n");
	    exit(4);
	}
    config[9] = (double)order;

#ifdef _WIN32
    printf("Configurating:\n");
#elif __linux__
    printf("\x1b[42;36mConfigurating:\x1b[0m\n");
#endif
    char * endptr;
    double conf_tmp;
    for (k = argc_least+1; k < argc; k++)
	{
	    errno = 0;
	    j = strtoul(argv[k], &endptr, 10);
	    if (errno != ERANGE && *endptr == '=')
		{							
		    endptr++;
		    errno = 0;
		    conf_tmp = strtod(endptr, &endptr);
		    if (errno != ERANGE && *endptr == '\0')
			{
			    config[j] = conf_tmp;
			    printf("%3d-th configuration: %g (ARGument)\n", j, conf_tmp);
			}
		    else
			{
			    printf("Configuration error in ARGument variable %d! ERROR after '='!\n", k);
			    exit(4);
			}
		}
	    else
		{
		    printf("Configuration error in ARGument variable %d! ERROR before '='!\n", k);
		    exit(4);
		}
	}
}
