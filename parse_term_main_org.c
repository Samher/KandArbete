#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include "parse_term_lib_org.c"
char *parse_vald_term(int, float, char *, char *, double);

int main(int npar, char *par[])
{
  char s[512], c_upp[3], c_low[3], level_low[87], level_upp[87], out[81], fname[23];
  char flag1, flag2, str1[4], str2[33]/*, wl[16]*/, *iret;
  float J_low, J_upp;
  double E_low, E_upp, wl, loggf;
  int species;
  FILE *fptr;

  iret=fgets(s, 511, stdin);
  iret=fgets(s, 511, stdin);
  
  //printf("par=%s\n", par[0]);
  //sprintf(fname, "parse_output_%c", );
  fptr = fopen("parse_output.txt","w");
  while(fgets(s, 511, stdin) != NULL)
  {
    //strncpy(wl, s, 16);
    wl=atof(s);
    J_low=atof(s+58);
    J_upp=atof(s+78);
    E_low=atof(s+44);
    E_upp=atof(s+66);
    species=atoi(s+30);
    loggf=atof(s+36);


    strncpy(level_low, s+129, 86); level_low[86]='\0';
    strncpy(level_upp, s+217, 86); level_upp[86]='\0';
    strncpy(c_low,s+127, 2); c_low[2]='\0';
    strncpy(c_upp,s+215, 2); c_upp[2]='\0';
    flag1=s[311];
/* Compare excitation energies to detect autoionisation lines */
    if(flag1 == 'A') flag1=(strncmp(s+44, s+64, 14)>0)?'A':' ';
    flag2=s[312];
//    printf("Level low%s\n", level_low);
//    printf("Coupling low=%s\n", c_low);
//    printf("Level upp=%s\n", level_upp);
//    printf("Coupling upp=%s\n", c_upp);
//    printf("Flag1=%c\n", flag1);

    strcpy(out,parse_vald_term(species, J_low, c_low, level_low, E_low));
    fprintf(fptr,"%s\n", out);
    

    strcpy(out,parse_vald_term(species, J_upp, c_upp, level_upp, E_upp));
    fprintf(fptr,"%s\n", out);

    switch(flag1)
    {
      case 'A': strcpy(str1, "A");
                break;
      case 'B': strcpy(str1, "E2");
                break;
      case 'C': strcpy(str1, "M1");
                break;
      case 'D': strcpy(str1, "M2");
                break;
      case 'E': strcpy(str1, "E3");
                break;
      case 'F': strcpy(str1, "M3");
                break;
      default:  str1[0]='\0';
                break;
    }

    switch(flag2)
    {
      case ' ': strcpy(str2, "none");
                break;
      case '0': strcpy(str2, "none");
                break;
      case '1': strcpy(str2, "Waals");
                break;
      case '2': strcpy(str2, "Stark");
                break;
      case '3': strcpy(str2, "Waals+Stark");
                break;
      case '4': strcpy(str2, "hfs");
                break;
      case '5': strcpy(str2, "Waals+hfs");
                break;
      case '6': strcpy(str2, "Stark+hfs");
                break;
      case '7': strcpy(str2, "Waals+Stark+hfs");
                break;
      default:  strcpy(str2, "none");
                break;
    }

    if(flag1 == 'A')
    {
      //printf("Autoionization, Extra_info: %s, %s\n", wl, str2);
      fprintf(fptr,"Autoionization,%.8lf,%.3lf,,%s\n", wl, loggf, str2);
    }
    else if(flag1 == ' ')
    {
      //printf("Allowed_transition: %s, E1, Extra_info: %s\n", wl, str2);
      fprintf(fptr,"Allowed_transition,%.8lf,%.3lf,E1,%s\n", wl, loggf, str2);
    }
    else
    {
      //printf("Forbidden_transition: %s, %s, Extra_info: %s\n", wl, str1, str2);
      fprintf(fptr,"Forbidden_transition,%.8lf,%.3lf,%s,%s\n", wl, loggf, str1, str2);
    }
    iret=fgets(s, 511, stdin);

  }
  fclose(fptr);
  return 0;
}
