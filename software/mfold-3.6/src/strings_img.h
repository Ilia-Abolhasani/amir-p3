#include <time.h>
#include <gd.h>
#include <gdfontg.h>
#include <gdfontl.h>
#include <gdfontmb.h>
#include <gdfonts.h>
#include <gdfontt.h>
#include <stdarg.h>

extern gdFontPtr gdFontTiny;
extern gdFontPtr gdFontSmall;
extern gdFontPtr gdFontMediumBold;
extern gdFontPtr gdFontLarge;
extern gdFontPtr gdFontGiant;

FILE *giffp;

gdImagePtr image;
gdImagePtr img_brush[15];

/* Variable format */
void gdImagePrintf(gdImagePtr im, gdFontPtr f, int x, int y, int color, 
		   char* s, ...) {
  char str[200];
  va_list ap;
  va_start(ap, s);
  vsnprintf(str, 200 , s, ap);
  va_end(ap);
  gdImageString(im, f, x, y, (unsigned char*) str, color);
}

void gdImageStringCenter(gdImagePtr image, gdFontPtr font, int x, int y, 
			 char* s, int color) {
  gdImageString(image, font, x - font->w * strlen(s)/2, y-font->h / 2, 
		(unsigned char*) s, color);
}

void gdImagePrintfCenter(gdImagePtr image, gdFontPtr font, int x, int y, 
			 int color, char* s, ...) {
  char str[200];
  va_list ap;
  va_start(ap, s);
  vsnprintf(str, 200 , s, ap);
  va_end(ap);
  gdImageString(image, font, x - font->w * strlen(str)/2, y - font->h / 2, 
		(unsigned char*) str, color);
}

/* Prototypes 
   void gdImagePrintf(gdImagePtr, gdFontPtr, int, int, int, char, ...);
   void gdImageStringCenter(gdImagePtr, gdFontPtr, int, int, char*, int);
   void gdImagePrintfCenter(gdImagePtr, gdFontPtr, int, int, int, char*, ...);
*/

void display_header_message_img(gdImagePtr image, int width, int annotation) {
  /*  char regc[4] = "(©)"; */
  char regc[4] = "(C)";
  char time_data[50] ;
  int length, headr, hor_pos, ver_pos;
  time_t now;
  now = time(NULL);
  strcpy(time_data, "Created ");
  strcat(time_data, ctime(&now));
  length = strlen(time_data);
  time_data[length-1] = '\0';
  hor_pos = width - 6*length; /* Treat as 6 point font */ 
  ver_pos = 1;
  /* Try not to make text overwrite image below */
  headr = gdImageColorAllocate(image, 0, 0, 50);
  gdImagePrintf(image, gdFontSmall, 6, 1, headr, "Output of %s %s", 
		PROGRAM_NAME, regc);
  gdImageString(image, gdFontSmall, 6, 15, 
		(unsigned char*) PACKAGE_STRING, headr);
  gdImagePrintf(image, gdFontSmall, hor_pos, ver_pos, headr, "%s", time_data); 
#if defined(NONE) && defined(P_NUM) && defined(SS_COUNT) && defined(PROB) 
  if (annotation != NONE) {
    if (annotation == P_NUM)
      gdImageString(image, gdFontSmall, 6, 29, 
		    (unsigned char*) "p-num annotation", headr);
    else if (annotation == SS_COUNT)
      gdImageString(image, gdFontSmall, 6, 29, 
		    (unsigned char*) "ss-count annotation", headr);
    else
      gdImageString(image, gdFontSmall, 6, 29, 
		    (unsigned char*) "probability annotation", headr);
  }
#endif
}
