void chop_suffix(char *full_string, char *mode ) {
  int i ;
  i = strlen(full_string);
  /* strip off ".gif", ".jpg" or ".png"  */
  if ( !strcmp(mode, ".img") ) {
    if (i > 4) {
      if ( !strcmp(full_string + i - 4, ".gif") || 
	   !strcmp(full_string + i - 4, ".jpg") || 
	   !strcmp(full_string + i - 4, ".png") )
	full_string[i - 4] = '\0';
    }
  } else {
    if (i > strlen(mode) ) {
      if ( !strcmp(full_string + i - strlen(mode), mode) )
	full_string[i - strlen(mode)] = '\0';
    }
  }
  return ;
}

void try_exit(int exit_choice) {
  if(exit_choice==0) {
    printf("Normal exit from %s.\n", PROGRAM_NAME);
    exit(0);
  } else {
    printf("Abnormal exit from %s.\n", PROGRAM_NAME);
    exit(exit_choice);
  }
}
