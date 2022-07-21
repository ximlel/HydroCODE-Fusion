#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>




/* This function examine whether a string
 * represents a real number.
 * Transform the string represents a
 * negtive number into a string represents
 * a positive one and return its' sign.
 * It returns 0 if the string do not
 * represents a real number.
 * After calling this function, there will
 * be only one 'e' in the string, and the
 * only position for '-' is behind 'e', and
 * there can be only one dot in the string
 * and the only position for it in before 'e'.
 */
int format_string(char * str)
{
  int i = 0, length = 0, j = 0;
  int sign = 1;
  int flag_dot = 0; //the number of dots in the string
                    //should be at most one
  int pos_dot = 0;
  int flag_e = 0;
  int pos_e = 0;

  //printf("%send\n", str);

  length = strlen(str);

  for(j = 0; j < length; ++j)
  {
    if((str[j] == 69) || (str[j] == 101))
    {
      str[j] = 101;
      flag_e += 1;
      pos_e = j;
    }
  }

  //there could not be more than one 'e' in one number
  if(flag_e > 1)
    return 0;
  if((flag_e) && (pos_e == 0))
    return 0;
  if((flag_e) && (pos_e == length-1))
    return 0;
  // a dot only could not be a number
  if((str[0] == 46) && (length == 1))
    return 0;
  // a '-' only could not be a number
  if(str[0] == 45)
  {
    if(length == 1)
      return 0;
    sign = -1;
  }

  // eliminate '-' from the string and return -1
  if(sign < 0)
  {
    for(i = 0; i < length; ++i)
      str[i] = str[i+1];
    length -= 1;
    pos_e -= 1;
    if(pos_e == 0)
      return 0;
  }

  if(flag_e)
  {
    for(i = 0; i < length; ++i)
    {
      if(str[i] == 45)
      {
      // after eliminate '-', the only possible position for '-'
      // is behind 'e'
        if((i-pos_e) != 1)
	  return 0;
	else if(i == length-1)
	  return 0;
      }
      //there could not be two dots in one number
      else if((str[i] == 46) && (flag_dot > 0))
        return 0;
      else if(str[i] == 46)
      {
        flag_dot += 1;
        pos_dot = i;
      }
    }
    if((flag_dot) && (pos_dot >= (pos_e-1)))
      return 0;
  }
  else
  {
    for(i = 0; i < length; ++i)
    {
      if(str[i] == 45)
        return 0;
      else if((str[i] == 46) && (flag_dot > 0))
        return 0;
      else if(str[i] == 46)
        flag_dot += 1;
    }
  }

  return sign;
}

/* this function trans form a string
 * consisting '1', '2', ..., and '.'
 * into the real number it represent
 */
double str2num(char * number)
{
  double result = 0.0, super_script = 0.0;
  int idx = 0, dot = -2;
  int i = 0, j = 0, power, k = 0;
  int length = 0;
  int pos_e = 0;
  char * after_e = number;
  int sign = 1;

  length = strlen(number);

  for(j = 0; j < length; ++j)
    if(number[j] == 101)
      pos_e = j;

  /*
  for(k = 0; k < length; ++k)
    printf("%c", number[k]);
  printf("\t");
  */

  if(pos_e)
  {
    after_e = number + pos_e + 1;
    number[pos_e] = 0;
    result = str2num(number);
    if(after_e[0] == 45)
    {
      sign = -1;
      after_e += 1;
    }
    super_script = str2num(after_e);
    result = result * pow(10.0, sign * super_script);
  }
  else
  {
    while(number[idx] != 0)
    {
      if(number[idx] == 46)
      {
        dot = idx - 1;
	idx = 0;
	break;
      }
      ++idx;
    }

    if(dot == -2)
      dot = idx - 1;

    for (i = 0; i <= dot; ++i)
      result += (double)(number[i] - 48)*pow(10, dot - i);

    dot += 1;
    double addon = 0.0;
    for (i = 1; i < length - dot; ++i)
      result += (double)(number[dot+i] - 48)*pow(0.1, i);
  }
  
  return result;
}


/* This function reads the initial data file
 * to generate the initial data.
 * It returns 0 if successfully read the file,
 * while returns the index of the wrong entry.
 */
int file_read(FILE * fp, double * U, int num)
{
  int idx = 0, j = 0; //j is a frequently used index for spatial variables
  char number[100];
  char ch;
  int sign = 1;

  while((ch = getc(fp)) != EOF)
  {
    //if(((ch <45) || ((ch > 57) && (ch != 69) && (ch != 101)) || (ch == 47)) && (idx))
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (idx))
    {
      number[idx] = 0;
      idx = 0;
      
//printf("%s\t",number);

      sign = format_string(number);
      if(!sign)
	return j+1;
      if(j == num)
	return j;

//printf("%d\t%s\n", sign, number);

      U[j] = sign * str2num(number);
      ++j;
    }
    else if((ch == 46) || (ch == 45) || ((ch > 47) && (ch < 58)) || (ch == 69) || (ch == 101))
      number[idx++] = ch;
  }

  return 0;
}
