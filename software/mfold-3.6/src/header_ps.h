#include <time.h>

#define TRUE    1
#define FALSE   0

void extra_ps_char(FILE* postfp) {
  /* for registered, copyright, trademark, delta, degree */
  fprintf(postfp,"\n/ISOLatin1Encoding where {pop save true}{false} ifelse");
  fprintf(postfp,"\n/reencodeISO {"); 
  fprintf(postfp,"\n   dup length dict begin");
  fprintf(postfp,"\n        {1 index /FID ne {def}{pop pop} ifelse} forall"); 
  fprintf(postfp,"\n         /Encoding ISOLatin1Encoding def");
  fprintf(postfp,"\n        currentdict"); 
  fprintf(postfp,"\n    end");
  fprintf(postfp,"\n} def"); 
  fprintf(postfp,"\n/findISO {");
  fprintf(postfp,"\n    dup /FontType known {"); 
  fprintf(postfp,"\n        dup /FontType get 3 ne {");
  fprintf(postfp,"\n              dup /CharStrings known {"); 
  fprintf(postfp,"\n                 dup /CharStrings get /Thorn known {");
  fprintf(postfp,"\n                     true"); 
  fprintf(postfp,"\n                }{ false } ifelse");
  fprintf(postfp,"\n            }{ false } ifelse");
  fprintf(postfp,"\n       }{ false } ifelse");
  fprintf(postfp,"\n     }{ false } ifelse"); 
  fprintf(postfp,"\n} def");
  fprintf(postfp,"\n");
  fprintf(postfp,"\n/delta { /Symbol findfont findISO { ");
  fprintf(postfp, "reencodeISO /Symbol-ISO exch definefont");
  fprintf(postfp,"\n                    sf scalefont setfont }");
  fprintf(postfp,"\n                   { sf scalefont setfont } ifelse");
  fprintf(postfp,"\n  (d) show } def");
  fprintf(postfp,"\n");
}

void display_header_message_ps(FILE *postfp, int portrait_flag) {
  char time_data[50];
  int length, cprt = 169;
  float hor_pos, ver_offset;
  time_t now;
  now = time(NULL);
  strcpy(time_data, "Created ");
  strcat(time_data, ctime(&now));
  length = strlen(time_data);
  time_data[length-1] = '\0' ;
    if(portrait_flag==TRUE) {
      ver_offset = 710.;
      hor_pos = 612 - 7.5*(float)length/2.0 - 15 ; /* 7.5 is font size */
    } else {
      ver_offset = 560;
      hor_pos = 792 - 7.5*(float)length/2.0 - 15 ;
    }
  extra_ps_char(postfp);
  /* Try not to make text overwrite image below */
  fprintf(postfp,"/sf 7.5 def\n");
  fprintf(postfp,"/Helvetica findfont findISO { reencodeISO /Symbol-ISO exch definefont \n");
  fprintf(postfp,"               sf scalefont setfont }\n");
  fprintf(postfp,"               { sf scalefont setfont } ifelse\n");
  fprintf(postfp,"15 %f moveto", ver_offset+66.);
  fprintf(postfp," (Output of %s (%c)) show\n" , PROGRAM_NAME, cprt);
  fprintf(postfp,"15 %f moveto", ver_offset+55);
  fprintf(postfp," (%s) show\n", PACKAGE_STRING);
  fprintf(postfp,"%f %f moveto", hor_pos, ver_offset+66);
  fprintf(postfp," (%s) show\n", time_data);
}
