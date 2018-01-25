

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <data.h>
#include <velocity_verlet.h>



/* main */
int main(int argc, char **argv)
{

    // char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp;

    /* pre-filled sys with selected values */
    int test1 = 3;
    double test2 = 39.948;
    char test3[BLEN] = "Test string";


    char inp1[BLEN];
    char inp2[BLEN];
    char inp3[BLEN];



    /* read restart */
    fp=fopen("test.dat","r");
    if(fp) {


     /* read input file */
     //if(get_a_line(stdin,line)) return 1;

     //integer case

	printf("Input int\n"); //atof
	if (get_a_line(fp,inp1)) {printf("error reading from file.\n"); return 1;}
        if (atoi(inp1)==test1) { printf("first input ok\n");}
	else  return 2;

     //floating case

	printf("Input float\n"); //atof
	if (get_a_line(fp,inp2)) {printf("error reading from file.\n"); return 1;}
	if(atof(inp2)==test2) { printf("second input ok\n"); }
	else  return 2;
     //string case
	
	printf("Input string\n"); //strcmp
	if (get_a_line(fp,inp3)) {printf("error reading from file.\n"); return 1;}
	if(!strcmp(inp3,test3)) { printf("third input ok\n"); }
	else  return 2;
	printf("Tests passed!\n");

	fclose(fp);


      }else{
		printf("No file...\n");
		return 1;
      }

	return 0;

      }


        
