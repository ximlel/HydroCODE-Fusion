/**
 * @file tools.h
 * @brief This file is the header file of several independent tool functions.
 * @details This header file declares functions in files 'sys_pro.c' and 'math_algo.c',
 */

#ifndef TOOLS_H
#define TOOLS_H

void DispPro(double pro, int step);
int CreateDir(const char* pPath);

int rinv(double a[], int n);
void Gauss_elimination(int n, double (*a)[n+1], double *x);

#endif
