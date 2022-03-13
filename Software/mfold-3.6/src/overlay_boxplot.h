
#define MAXIMUM_HELICES 400000 

#define MAXIMUM_CT_FILES 1 /* Limit on number of ct files that can be
			      read. Do not change it for this application */

#define MAXIMUM_SIZE 30000 /* maximum sequence length */
#define TRUE    1
#define FALSE   0

#define MAIN_POST    1 /* create postcript for main or zoom */
#define ZOOM_POST    2 /* used in general_post, create_post */

struct helix {
  int row;
  int column;
  int diagonal;
  int color;
  int energy;
  int length;
};

struct ct_helix {
  int row;
  int column;
  int diagonal;
  int color;         
  int length;
  int plotlen;
  int energy;
};

/* length is length of ct helix at fixed energy within plot file
 * plotlen is the length of the helix in the plot file that has this
 * energy stored as energy here 
 */

float g_close_optimal; /* defines what is shown as yellow */
int energy_read(int,int);
int chain_length_read(int,int);
int find_color(int);

/*__________________________tick mark functions ___________________
 *
 * Determine distance between tick marks */
int step_fun(int dis_l,int dis_r) {
  int step;
  step=(dis_r-dis_l+1)/8; /* make 8 steps  (9 tick marks total) */
  if(step<40) {
    if(step<3)
      step=2;
    else if (step>5) {
      step = (step+5)/10; /* make a large step a multiple of ten */
      step = step*10;
    } else if (step==4) /* make steps of 4 look like 5 */
      step=5;
  } else if(step>90) {
    step=(step+50)/100; /* round to nearest 100 */
    step=step*100;
  } else {
    step = (step+25)/50; /* round to nearest 50 */
    step = step*50;
  }
  return step;
}

/* make the tick marks start at a nice place.
 * ends are treated separately and always drawn */
int start_fun(int dis_l,int step) {
  int start;
  start = dis_l + step;
  if(step>90) {
    start = (start+60)/100; /* start at a multiple of 100 */
    start = start*100;
  } else if(step>=40) {
    start=(start+32)/50; /* start at a multiple of 50 */
    start=start*50;
  } else if(step>5) {
    start=(start+5)/10; /* start at a multiple of ten */
    start=start*10;
  }
  return start;
}

/* Conversion Routines to display integers and floats */
char g_str[30];
char *num_string_int(int x) { /* converts integer to string */
  sprintf(g_str,"%d",x);
  return g_str;
}

char *num_string_float(float x) { /* convert float to string */
  sprintf(g_str,"%.3f",x);
  return g_str;
}

char *num_prob_string(int x) { /* divides by 1000 and multiplies by -1 */
  float x_float;
  x_float = (float) x;
  x_float = x_float/1000.;
  x_float = -x_float;
  /* round to thousandth position */
  sprintf(g_str,"%.3f",x_float);
  return g_str;
} 

/* divides number by 10 and converts to string */
char *num_energy_string(int x) {
  float x_float;
  x_float=(float)x;
  x_float=x_float/10.;
  /* round to tenth position */
  sprintf(g_str,"%7.1f",x_float);
  return g_str;
}

char *num_prob_string_float(float x) {
  x = -x/1000.;
  /* round to thousandth position */
  sprintf(g_str,"%4.3f",x);
  return g_str;
} 

char *num_energy_string_float(float x) {
  x = x/10.;
  /* round to tenth position */
  sprintf(g_str,"%7.1f",x);
  return g_str;
}

char *num_string_fancy_int(int x,int prob_flag) {
  /* converts integer to string , handles prob or energies*/
  float float_x;
  float_x=(float)x;
  if(prob_flag)
    return num_prob_string_float(float_x);
  else
    return num_energy_string_float(float_x);
}

/* converts float to string */
char *num_string_fancy_float(float x,int prob_flag)  {
  if(prob_flag)
    return num_prob_string_float(x);
  else
    return num_energy_string_float(x);
}

int set_ct_color(int energy,int optimal_energy,int worst_energy,
                 float close_optimal_energy) {
  /* return 0 for optimal, 1 for close, 2 for not optimal */
  if(energy==optimal_energy) return 0;
  if(energy==INT_MIN) return 2;
  if(worst_energy==optimal_energy) return 0;
  if(((float)(energy-optimal_energy))/
      ((float)(worst_energy-optimal_energy))<close_optimal_energy)
    return 1;
  else 
    return 2;
}
