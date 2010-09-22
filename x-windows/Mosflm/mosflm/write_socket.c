/* 1 Make the necessary includes and set up the variables: */

#include <sys/types.h>
#include <stdio.h>
#if (defined _MSC_VER) || (defined __WIN32)
#include <winsock.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

int write_socket_length_(
  int *fd,
  int *length,
  char *line
)
{
#if (defined _MSC_VER) || (defined __WIN32)
  send(*fd, line, *length,(int) NULL);
  send(*fd, "\n", 1,(int) NULL);
#else
  write(*fd, line, *length);
  write(*fd, "\n", 1);
#endif

  return 0;
}

int write_socket_section_(
  int *fd,
  int *length,
  char *line
)
{
#if (defined _MSC_VER) || (defined __WIN32)
  send(*fd, line, *length,(int) NULL);
#else
  write(*fd, line, *length);
#endif
  return 0;
}
