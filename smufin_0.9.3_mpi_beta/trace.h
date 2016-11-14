#ifndef __TRACE_H
#define __TRACE_H

#ifdef HAVE_LOG4C_H
#include <log4c.h>
#else

#define LOG4C_PRIORITY_TRACE 10
#define LOG4C_PRIORITY_DEBUG 20

#endif

#ifndef LOGGERNAME
#define LOGGERNAME "root"
#endif

#ifdef LOGGING_ENABLED

void TraceMsg(char *loggerName, int level, char *format, ...);
void TraceEnter(char *loggername, const char *function, char *format, ...);
void TraceLeave(char *loggername, const char *function, char *format, ...);
void TraceShowIntString(char *loggerName, int level, const int *str, int start, int len);

#define LOG(level,...) TraceMsg(LOGGERNAME, level, __VA_ARGS__);
#define SHOWINTSTRING(level,str,start,len) TraceShowIntString(LOGGERNAME, level, str, start, len);
#define ENTER(...) TraceEnter(LOGGERNAME,__PRETTY_FUNCTION__,__VA_ARGS__);
#define LEAVE(...) TraceLeave(LOGGERNAME,__PRETTY_FUNCTION__,__VA_ARGS__);

#else

#define LOG(level,...)
#define SHOWINTSTRING(level,str,start,len)
#define ENTER(...)
#define LEAVE(...)

#endif

#endif
