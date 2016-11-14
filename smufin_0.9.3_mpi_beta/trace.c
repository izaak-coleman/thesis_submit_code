/*
 * trace.c
 */

#include <string.h>

#ifdef HAVE_LOG4C
#include <log4c.h>
#endif

#include "trace.h"


#define MIN(a,b) ((a)<(b)?(a):(b))
#define LOG_PRI 0

#ifdef LOGGING
static int trace_initialized;
#endif

void TraceMsg(char *loggerName, int level, char *format, ...)
{
#ifdef LOGGING
	va_list args;
	log4c_category_t *cat;
	
	cat = log4c_category_get(loggerName);

	if (!trace_initialized)
	{
		log4c_category_set_appender(cat,log4c_appender_get("stdout"));
		log4c_category_set_priority(cat,LOG_PRI);
		trace_initialized = 1;
	}

	va_start(args,format);
	log4c_category_vlog(cat, level, format, args);
	va_end(args);
#endif
}

void TraceEnter(char *loggerName, const char *function, char *format, ...)
{
#ifdef LOGGING
	log4c_category_t *cat;
	
	char enterBuf[1024];
	char newformat[200];

	va_list list;
	va_start(list, format);

	strcpy(newformat,"> Enter ");
	strlcat(newformat,function,sizeof(newformat));
	strlcat(newformat,"(",sizeof(newformat));
	strlcat(newformat,format,sizeof(newformat));
	strlcat(newformat,")",sizeof(newformat));

	vsnprintf(enterBuf,sizeof(enterBuf),newformat,list);

	va_end(list);

	cat = log4c_category_get(loggerName);
	if (!trace_initialized)
	{
		log4c_category_set_appender(cat,log4c_appender_get("stdout"));
		log4c_category_set_priority(cat,LOG_PRI);
		trace_initialized = 1;
	}

	log4c_category_log(cat, LOG4C_PRIORITY_TRACE, enterBuf);
#endif
}

void TraceLeave(char *loggerName, const char *function, char *format, ...)
{
#ifdef LOGGING
	log4c_category_t *cat;
	
	char enterBuf[1024];
	char newformat[200];

	va_list list;
	va_start(list, format);

	strcpy(newformat,"< Leave ");
	strlcat(newformat,function,sizeof(newformat));
	strlcat(newformat,"(",sizeof(newformat));
	strlcat(newformat,format,sizeof(newformat));
	strlcat(newformat,")",sizeof(newformat));

	vsnprintf(enterBuf,sizeof(enterBuf),newformat,list);

	va_end(list);

	cat = log4c_category_get(loggerName);
	if (!trace_initialized)
	{
		log4c_category_set_appender(cat,log4c_appender_get("stdout"));
		log4c_category_set_priority(cat,LOG_PRI);
		trace_initialized = 1;
	}

	log4c_category_log(cat, LOG4C_PRIORITY_TRACE, enterBuf);
#endif
}

void TraceShowIntString(char *loggerName, int level, const int *str, int start, int len)
{
#ifdef LOGGING
	int i;
	char buf[260];
	log4c_category_t *cat;

	for (i=start;i<start+MIN(sizeof(buf)-4,len);i++)
		buf[i] = str[i];
	buf[i++] = '.';
	buf[i++] = '.';
	buf[i++] = '.';
	buf[i] = 0;

	cat = log4c_category_get(loggerName);
	if (!trace_initialized)
	{
		log4c_category_set_appender(cat,log4c_appender_get("stdout"));
		log4c_category_set_priority(cat,LOG_PRI);
		trace_initialized = 1;
	}

	log4c_category_log(cat, level, buf);
#endif
}
#ifdef TEST

int main(void)
{
	LOG(LOG4C_PRIORITY_FATAL,"exSuffixArray_create()");
}

#endif
