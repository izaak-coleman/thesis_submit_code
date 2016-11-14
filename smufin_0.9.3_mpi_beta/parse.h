#ifndef SUFFIX__PARSE_H
#define SUFFIX__PARSE_H

char *skip_spaces(char *line);
char *skip_simple_string(char *line);
char *parse_simple_string(char *line, char **dest_ptr);
char *parse_double(char *line, double *val_ptr);

#endif /*SUFFIX__PARSE_H*/
