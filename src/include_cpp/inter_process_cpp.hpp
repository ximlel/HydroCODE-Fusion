/**
 * @file inter_process_cpp.hpp
 * @brief This file is the c++ header file of intermediate processes of finite volume scheme.
 * @details This header file declares functions in the folder 'inter_process_cpp'.
 */

#ifndef _INTER_PROCESS_CPP_HPP_
#define _INTER_PROCESS_CPP_HPP_


#ifdef __cplusplus
extern "C" {
#endif
///////////////////////////////////
// VIPLimiter.cpp
///////////////////////////////////
double useVIPLimiter(const int neigh_cell_num, const double Vave[][2], double* V0, double* Vp);
#ifdef __cplusplus
}
#endif

#endif
