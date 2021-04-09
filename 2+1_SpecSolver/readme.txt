1) 'make' should compile the code
2) 'make clean' should clean executable and obj files

OBS: In debug.c there is a function called "pause()"
if I compile in MAC, I need to use fpurge(stdin);
if I compile in LINUX, I need to use __fpurge(stdin);

(if someone can code this if/else condition according to OS, it would be great....)