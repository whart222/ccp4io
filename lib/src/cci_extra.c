#include <stdio.h>
#include <string.h>

#if defined(GFORTRAN)
#include <time.h>
#endif

#if defined(CCP4LIB_CCI_IFORT)
#include <errno.h>
#endif

#if defined(GFORTRAN) /* for solve_resolve */

void time_(char* result, int len_result)
{
  struct tm *lt;
  time_t tim;
  if (len_result != 8) {
    fprintf(stderr, "ERROR: time_() IMPLEMENTATION LIMITATION:"
                    " len(argument) must be 8 exactly.\n");
    exit(1);
  }
  tim = time(NULL);
  lt = localtime(&tim);
  sprintf(result, "%02d:%02d:%02d", lt->tm_hour, lt->tm_min, lt->tm_sec);
}

void date_(char* result, int len_result)
{
  static const char* months[] = {
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  };
  struct tm *lt;
  time_t tim;
  if (len_result != 9) {
    fprintf(stderr, "ERROR: date_() IMPLEMENTATION LIMITATION:"
                    " len(argument) must be 9 exactly.\n");
    exit(1);
  }
  tim = time(NULL);
  lt = localtime(&tim);
  sprintf(result, "%02d-%s-%02d",
    lt->tm_mday, months[lt->tm_mon], lt->tm_year % 100);
}

int isatty_(int* lunit)
{
  fprintf(stderr, "WARNING: isatty_(%d) NOT SUPPORTED.\n", lunit);
  return 0;
}

#endif /* GFORTRAN */

#if defined(CCP4LIB_CCI_IFORT)

void gerror_(char* result, int len_result)
{
  int i;
  if (errno == 0) {
    for (i=0;i<len_result;i++) result[i] = ' ';
  }
  else {
    strncpy(result, strerror(errno), len_result);
    for (i=strlen(result);i<len_result;i++) result[i] = ' ';
  }
}

#endif /* CCP4LIB_CCI_IFORT */
