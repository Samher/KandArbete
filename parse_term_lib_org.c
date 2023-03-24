#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

char *strtrim(char *s)
{
  int i,j,l;

  l=strlen(s);
  for(i=0;i<l;i++)
  {
    if(!isblank(s[i])) break;
  }
  for(j=l-1;j>=i;j--)
  {
    if(!isblank(s[j])) break;
  }
  if(i>0) strncpy(s,s+i,j-i+i);
  s[j-i+1]='\0';
  return s;
}

char *strcompr(char *s)
{
  int i,j,l;

  l=strlen(s);
  for(i=j=0;i<l;i++)
  {
    if(!isblank(s[i])) s[j++]=s[i];
  }
  s[j]='\0';
  return s;
}

char *strcompr1(char *s)
{
  int i,j,l,blank;

  blank=0;
  l=strlen(s);
  for(i=j=0;i<l;i++)
  {
    if(!isblank(s[i]))
    {
      s[j++]=s[i];
      blank=0;
    }
    else if(isblank(s[i]) && !blank)
    {
      s[j++]=s[i];
      blank=1;
    }
  }
  s[j]='\0';
  return s;
}

int a2i_half(char *s)
{
  return (strchr(s, '/')==NULL)?atoi(s):atoi(s)/2.;
}

float a2f_half(char *s)
{
  return (strchr(s, '/')==NULL)?atof(s):atof(s)/2.;
}

double a2d_half(char *s)
{
  return (strchr(s, '/')==NULL)?atof(s):atof(s)/2.;
}

double max_test( double param1, double param2)
{
   if (param1 > param2)
     {
       return param1;  // Notice: that param1 is of type double and the return
                       //         type is also of type double
     }
   else
     {
       return param2;
     }
 }

char *parse_vald_term(int species, float J, char coupling[], char level_name[], double energy)
{
  int l;
  char c, *j1, *j2, *ss1;
  char static out[81];
  char term[87], level[87];
  int parity;
  
  strncpy(level, level_name, 87);
  strtrim(level);
  j2=strcompr1(level);
/*  j2=level; */
  j1=strstr(level,"\\ ");
  c=level[strlen(level)-1];
  parity=2;  // Even parity = 2, odd parity = 1
  if(c=='*') parity=1;

  ss1=strstr(level_name, "n=");
/*  if(ss1 != NULL) printf("term:%s\n",level_name); */
  if(species<=3 && ss1 != NULL) /* Hydrogen and Helium levels above 1 */
  {
    //sprintf(out,"Unknown coupling:J,n,parity:%4.1f,%d,%d",J,atoi(ss1+2),parity);
    sprintf(out,"%s,%lf,UC,%4.1f,%d,%d",level,energy,J,atoi(ss1+2),parity);
    return out;
  }

  while(j1!=NULL)
  {
    j2=j1+2;
    j1=strstr(j2,"\\ ");
  }
  j1=strchr(j2,' ');
  if(j1==NULL)
  {
    //sprintf(out,"Unknown coupling:J,parity:%4.1f,%d",J,parity);
    sprintf(out,"%s,%lf,UC,%4.1f,%d",level,energy,J,parity);
    return out;
  }
  strcpy(term,j1);
  strtrim(term);
  *j1='\0';



/* Level description in VALD:

Each level is characterized by energy E in cm-1, total angular
momentum J, electronic configuration and term designation (name).

Electronic configuration is separated from term name by one or
more spaces and both are packed in to a fixed length (46 bytes)
field. If a space occurs inside an electronic configuration,
it is marked by a leading back slash "\ " (see examples below).

Odd parity levels have * in term names. Other information
(quantum numbers) can be extracted from term names and electronic
configurations.

In VALD we consider 3 types of coupling schemes
with the notations adopted by the NIST (see Atomic Spectroscopy,
part 9, Notations for different coupling schemes)

I. LS coupling

In LS coupling term name contains one block letter other than A
and B. This block letter gives the value of the orbital quantum
number L (S P D F G H I K L M  N .... =
          0 1 2 3 4 5 6 7 8 9 10 ....)
A number in front of the block letter is multiplicity and it is
equal to 2S+1, where S is spin. Any small letters in front of the
multiplicity are ignored. Any numbers or block letters A, B after
block letter containing L information are seniority indices and
can be ignored when extracting quantum numbers.

Examples:

Ion     E          J         Electronic configuration        Term name
Tm I   8771.240  2.5   LS  4f13.(2F*).6s2                        2F*
Fe I  60757.592  4.0   LS  3p6.3d6.(3H).4s.4p.(1P*)             t3H*
Fe I  57070.167  6.0   LS  3p6.3d6.(1I).4s.4p.(3P*)             x3I*
Fe I  26874.548  5.0   LS  3p6.3d6.(5D).4s.4p.(3P*)             z5F*
Fe I  29313.006  6.0   LS  3p6.3d6.4s2                          a1I
Fe I  26623.733  2.0   LS  3p6.3d7.(a2D).4s                     a3D
Fe II 20516.960  2.5   LS  3p6.3d7                              a2D2

The first four levels have odd parity. For the last level number
"2" after "D" is seniority. All these levels have the following
sets of quantum numbers:
      J    L     S
1.   2.5   3    0.5
2.   4.0   5    1.0
3.   6.0   6    1.0
4.   5.0   3    2.0
5.   6.0   6    0.0
6.   2.0   2    1.0
7.   2.5   2    0.5

L, S vectors are coupled to give J = from L-S to L+S.
*/

  if(!strncmp(coupling,"LS",2)) /* LS coupling, parse quantum numbers */
  {
    int L;              /* Orbital quantum number                   */
    float S;            /* Spin                                     */
    int multiplicity;
    char seniority;
    int decimal, ll;

    l=strlen(term)-1;
    c=term[l];

    parity=2;
    seniority=' ';
    L=0;
    if(c=='*') {parity=1; l--; c=term[l];} /* Parse and skip parity */
    if(c=='?') {l--; c=term[l];} /* Skip leftovers */
//    if(c>='A' && c<='z') {l--; c=term[l];} /* Skip leftovers */
    if(c=='*') {parity=1; l--; c=term[l];} /* Parse and skip parity if it shows up again */
    if(c=='X') {l--; c=term[l];} /* Skip strange X */
    if(c=='a') {l--; c=term[l];} /* Skip strange a */
    if(c=='b') {l--; c=term[l];} /* Skip strange b */
    if(c=='c') {l--; c=term[l];} /* Skip strange c */
    if(c=='+') {l--; c=term[l];} /* Skip strange + */

    if((c<='9' && c>='0') || c=='A' || c=='B') {seniority=c; l--; c=term[l];}

//  printf("term: %s\n", level_name);
//  printf("term: %s\n seniority:%c, parity:%d, l=%d\n", term, seniority, parity, l);

    ll=l-1;
    multiplicity=0;
    decimal=1;
    while(ll>=0 && isdigit(term[ll]))
    {
      multiplicity+=((int)(term[ll]-'0'))*decimal;
      ll--;
      decimal*=10;
    }
    S=(multiplicity-1)/2.;

    if(c=='S')      L= 0;
    else if(c=='P') L= 1;
    else if(c=='D') L= 2;
    else if(c=='F') L= 3;
    else if(c=='G') L= 4;
    else if(c=='H') L= 5;
    else if(c=='I') L= 6;
    else if(c=='K') L= 7;
    else if(c=='L') L= 8;
    else if(c=='M') L= 9;
    else if(c=='N') L=10;

    if(seniority==' ')
    {
/*    sprintf(out,"LS,J:%4.1f,L:%d,S:%4.1f,parity:%d",J,L,S,parity); */
      //sprintf(out,"Label:%s,E:%f,LS:J,L,S,parity:%4.1f,%d,%4.1f,%d",level,E,J,L,S,parity);
      sprintf(out,"%s,%lf,LS,%4.1f,%d,%4.1f,%d",level,energy,J,L,S,parity);
    }
    else
    {
/*    sprintf(out,"LS,J:%4.1f,L:%d,S:%4.1f,parity:%d,seniority:%c",J,L,S,parity,seniority); */
      //sprintf(out,"Label:%s,E:%f,LS:J,L,S,parity,seniority:%4.1f,%d,%4.1f,%d,%c",level,E,J,L,S,parity,seniority);
      sprintf(out,"%s,%lf,LS,%4.1f,%d,%4.1f,%d,%c",level,energy,J,L,S,parity,seniority);
    }

    if(S<0)
    {
      printf("Negative S=%g, mutiplicity=%d decimal=%d term=%s\n",S,multiplicity,decimal,term);  // How to handle this?
      printf("%s\n", out);
    }

    strcompr(out);
    return out;
  }

/*
II. JJ coupling

In JJ coupling scheme term name contains two J values in (J1,J2),
separated by comma. It is usual that J1 and J2 are also given
in electronic configuration in <>. In this coupling scheme each
level has the following quantum numbers: J, J1 and J2

Examples:

Ion    E          J         Electronic configuration       Term name
Tm I 37858.620  7.500 JJ  4f12.(3H<6>).5d.6s.6p.(4D*<7/2>) (6,7/2)*

J=7.5, J1=6, J2=3.5

J1, J2 are coupled to give J = abs(J1-J2), ..., J1+J2.
Number of levels is 2*min(J1,J2)+1.

Example: 3s2.3p5.(2P*<1/2>).5s (1/2,1/2)*
*/

  if(!strncmp(coupling,"JJ",2)) /* JJ coupling, parse quantum numbers */
  {
    char t[60];
    float J1;
    float J2;

    l=strlen(term)-1;
    c=term[l];
    parity=2;
    if(c=='*') {parity=1; l--; c=term[l];}
    J1=0.; J2=0.;
    j1=strchr(term,'('); j2=strchr(term,',');
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,(int)(j2-j1)-1);
      t[((int)(j2-j1))-1]='\0';
      J1=a2f_half(t);
    }
    j1=strchr(term,','); j2=strchr(term,')');
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,(int)(j2-j1)-1);
      t[((int)(j2-j1))-1]='\0';
      J2=a2f_half(t);
    }
/*  sprintf(out,"JJ,J:%4.1f,J1:%4.1f,J2:%4.1f,parity:%d",J,J1,J2,parity); */
    //sprintf(out,"Label:%s,E:%f,JJ:J,J1,J2,parity:%4.1f,%4.1f,%4.1f,%d",level,E,J,J1,J2,parity);
    sprintf(out,"%s,%lf,JJ,%4.1f,%4.1f,%4.1f,%d",level,energy,J,J1,J2,parity);
/*    printf("JJ:J,J1,J2,parity:%4.1f,%4.1f,%4.1f,%d",J,J1,J2,parity); */
    strcompr(out);
    return out;
  }

/*
III. JK coupling

In JK coupling scheme each level is characterized by the
following set of quantum numbers: J, Jc, S2 and K. The latter two
can be extracted from term name, where number in [] gives K, and
number in front of [] is 2*S2+1. Jc is J-value of the parent term
and is usually given in <> inside electronic configuration.

Examples:

Ion      E       J          Electronic configuration      Term name
Fe I 61541.170  5.0   JK  3p6.3d6.(5D).4s.\ (6D<3/2>).6g   2[11/2]
Tm I 39322.014  4.5   JK  4f13.(2F*<7/2>).6s.6d.(3D)        3[7/2]*

K and S2 are coupled to give J: J=K+/-S2, the multiplicity of [K]
term arises from the spin of the external electron(s). In first
example we have one external electron, 6g, therefore S2=0.5, and
each term produces pairs of levels with J=K+/-0.5 and we have two
levels with J=5.0 and 6.0. Here Jc=1.5 = 3/2. In the second
example S2=1 (spin of 3D), Jc=3.5, K=3.5, and this term produces
3 levels with J=4.5, 3.5, 2.5

A space inside electronic configuration is noted as "\ ", and it
emphasize the order of coupling: 4s electron first is coupled to
5D grandparent term to produce 6D<3/2>, and then 6g electron is
coupled to 6D<3/2> parent to give a resulting term 2[11/2].
*/

  if(!strncmp(coupling,"JK",2)) /* JK coupling, parse quantum numbers */
  {
    char t[60];
    float K, Jc, S2;

    l=strlen(term)-1;
    c=term[l];
    parity=2;
    K=0.;
    S2=0.;
    Jc=0.;
    if(c=='*') {parity=1; l--; c=term[l];}
    j1=strchr(term,'('); j2=strchr(term,')');
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,((int)(j2-j1))-1);
      t[((int)(j2-j1))-1]='\0';
      S2=a2f_half(t);
    }
    else
    {
      S2=atof(term);
      S2=(S2-1.)/2.;
    }
    j1=strchr(term,'['); j2=strchr(term,']');
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,((int)(j2-j1))-1);
      t[((int)(j2-j1))-1]='\0';
      K=a2f_half(t);
    }
    j1=strchr(level,'<'); j2=strchr(level,'>');
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,((int)(j2-j1))-1);
      t[((int)(j2-j1))-1]='\0';
      Jc=a2f_half(t);
    }
/*  sprintf(out,"JK,J:%4.1f,K:%4.1f,Jc:%4.1f,S2:%4.1f,parity:%d",J,K,Jc,S2,parity); */
    if(Jc>0.01)
    {
      /*sprintf(out,"Label:%s,E:%f,JK:J,K,Jc,S2,parity:%4.1f,%4.1f,%4.1f,%4.1f,%d",
                      level,E,J,K,Jc,S2,parity);*/
      sprintf(out,"%s,%lf,JK,%4.1f,%4.1f,%4.1f,%4.1f,%d",
                      level,energy,J,K,Jc,S2,parity);
    }
    else
    {
      /*sprintf(out,"Label:%s,E:%f,JK:J,K,S2,parity:%4.1f,%4.1f,%4.1f,%d",
                      level,E,J,K,S2,parity);*/
      sprintf(out,"%s,%lf,JK,%4.1f,%4.1f,%4.1f,%d",
                      level,energy,J,K,S2,parity);
    }
    strcompr(out);
    return out;
  }

/*
IV. LK coupling

In LK coupling scheme each level is characterized by the following set
of quantum numbers: J, L, S2 and K. S2 and K can be extracted from term
name, where the number in square brackets gives K, and number in front
of the [] is 2*S2+1. L is extracted from the block letter in the end
of electronic configuration, separated from the configuration by an
escaped space "\ ".

Example:

Ion          E       J    Coupling   Electronic configuration   Term name
P II    132132.460  3.0     LK          3s2.3p.(2P*).4f.\ D       2[5/2]

Letter "D" gives us L=2 (sum of orbital angular momenta of the core and
spin of the external electron), while S2 and K are extracted from term
name as in JK coupling. Here vector K is a sum of L and S1 - spin of
the core (S1=1/2 in example).
*/
  if(!strncmp(coupling,"LK",2)) /* LK coupling, parse quantum numbers */
  {
    char t[60], c;
    float K, S2;
    int L;

    l=strlen(term)-1;
    c=term[l];
    parity=2;
    K=0.;
    S2=0.;
    L=-1;

    if(c=='*') {parity=1; l--; c=term[l];}
    j1=strchr(term,'('); j2=strchr(term,')'); /* Parsing term */
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,((int)(j2-j1))-1);
      t[((int)(j2-j1))-1]='\0';
      S2=a2f_half(t);
    }
    else
    {
      S2=atof(term);
      S2=(S2-1.)/2.;
    }
    j1=strchr(term,'['); j2=strchr(term,']');
    if(j1!=NULL&&j2>j1)
    {
      strncpy(t,j1+1,((int)(j2-j1))-1);
      t[((int)(j2-j1))-1]='\0';
      K=a2f_half(t);
    }

    j1=strstr(level,"\\ ");                   /* Parsing electronic config. */
    if(j1!=NULL)
    {
      c=*(j1+2);
      if(c=='S')      L= 0;
      else if(c=='P') L= 1;
      else if(c=='D') L= 2;
      else if(c=='F') L= 3;
      else if(c=='G') L= 4;
      else if(c=='H') L= 5;
      else if(c=='I') L= 6;
      else if(c=='K') L= 7;
      else if(c=='L') L= 8;
      else if(c=='M') L= 9;
      else if(c=='N') L=10;
    }

    if(L>=0)
    {
      /*sprintf(out,"Label:%s,E:%f,LK:J,K,L,S2,parity:%4.1f,%4.1f, %d,%4.1f,%d",
                      level,E,J,K,L,S2,parity);*/
      sprintf(out,"%s,%lf,LK,%4.1f,%4.1f, %d,%4.1f,%d",
                      level,energy,J,K,L,S2,parity);
    }
    else
    {
      /*sprintf(out,"Label:%s,E:%f,LK:J,K,S2,parity:%4.1f,%4.1f,%4.1f,%d",
                      level,E,J,K,S2,parity);*/
      sprintf(out,"%s,%lf,LK,%4.1f,%4.1f,%4.1f,%d",
                      level,energy,J,K,S2,parity);                
    }
    strcompr(out);
    return out;
  }

/*
V. Unknown coupling.

In this case there is no term name. Sometimes electronic configuration
exists. In this case Term Type = 'Unknown', and energy level is characterized
by energy value, J and parity.

****************************************************************************************************************
In VALD all J for individual levels are written in the form 1.0,
1.5, etc, while J-numbers inside configuration and in term names
(J1,J2,K) are usually written as 1, 3/2, etc, although sometimes
1, 1.5, ... may be used.
*/
  //sprintf(out,"Label:%s,E:%f,Unknown coupling:J,parity:%4.1f,%d",level,E,J,parity);
  sprintf(out,"%s,%lf,UC,%4.1f,%d",level,energy,J,parity);
  return out;
}
