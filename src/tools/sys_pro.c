/**
 * @file  math_algo.c
 * @brief There are some system processing programs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

/*
 * To realize cross-platform programming.
 * MKDIR:  Create a subdirectory.
 * ACCESS: Determine access permissions for files or folders.
 */
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#include <io.h>
/*
 * m=0: Test for existence.
 * m=2: Test for write permission.
 * m=4: Test for read permission.
 */
#define ACCESS(a,m) _access((a),(m))
#define MKDIR(a)    _mkdir((a))  // Create a subdirectory.
#elif __linux__
#include <sys/stat.h>
#include <unistd.h>
#define ACCESS(a,m) access((a),(m))
#define MKDIR(a)    mkdir((a),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif


void DispPro(double pro, int step)
{
        int j;
        for (j = 0; j < 77; j++)
                putchar('\b'); // 将当前行全部清空，用以显示最新的进度条状态
        for (j = 0; j < lround(pro/2); j++)
                putchar('+'); // 打印进度条上已经完成的部分，用‘+’表示  
        for (j = 1; j <= 50-lround(pro/2); j++)
                putchar('-'); // 打印进度条上还有多少没有完成的  
        fprintf(stdout, "  %6.2f%%   STEP=%-8d", pro, step);  
        fflush(stdout);
}

int CreateDir(const char * pPath)
{
	if(-1 != ACCESS(pPath,2))
		return -1;

	const char* pCur = pPath;

	char tmpPath[FILENAME_MAX];
	memset(tmpPath,0,sizeof(tmpPath));
    
	int pos = 0;
	while(*pCur++!='\0')
		{
			tmpPath[pos++] = *(pCur-1);

			if(*pCur=='/' || *pCur=='\0')
				{
					if(0!=ACCESS(tmpPath,0) && strlen(tmpPath)>0)
						{
							MKDIR(tmpPath);
						}
				}
		}
	if(!ACCESS(pPath,2))
		return 0;
	else
		return 1;
}
