/* temporary wrapper just to read Pilatus CBF image data */

#include <include/cbf.h>
#include <examples/img.h>
#include <include/cbf_simple.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#undef cbf_failnez
#define cbf_failnez(x) \
  {int err;	       \
    err = (x);				      \
    if (err) {					      \
    }							\
  }

int pilatus_c_(int *ierr, int *modeop, int *dim_1, int *dim_2, short image[64000000], char *argv) {

  FILE *in, *out;

  clock_t a,b;

  img_handle img, cbf_img;

  cbf_handle cbf;
  cbf_detector this_detector;
  cbf_goniometer this_goniometer;
  /* temporary from cbfwrap.c - pilatus doesn't have these, so we'll 
delete them later */
  int cbf_int[16];
  long cbf_int4[16];
  double cbf_double[16];
  char cbf_char[16][24];
  /* temporary from cbfwrap.c - pilatus doesn't have these, so we'll 
delete them later */
  int id, index, iindex;

  int i,ii,jjj,colrow, first,second,dirsta[2],dirend[2],dirinc[2],tmp_pixel;

  unsigned int column, row;

  size_t nelem_read;

  double pixel_size, gain, wavelength, distance, phi_start, phi_range;
  double beam_centre[4];
  int overload, dimension[2], precedence ;

  const char *detector,*oscillation_axis;

  char *detector_char;

  char detector_id [64];

  char *asl_index[2];

  const char *direction [2], *array_id, *radiation;
  
  int maxpix, minpix, minx, miny, maxx, maxy ;
  
  /* the following declarations are for transferring data between 
     this C routine and the FORTRAN calling routine */
  minpix = 999999;
  maxpix = 0;
  /* start time for file read */
  a = clock ();
                      
  direction[0] = "null";
  direction[1] = "null";
  array_id = "null";
  /* Create the cbf handle for the image file*/

  cbf_failnez (cbf_make_handle (&cbf))


    /* Read the file */

  in = fopen (argv, "rb");
  if (!in)
  {
    fprintf (stderr, " Couldn't open the CBF file %s\n", argv);
    *ierr = 1;
    return (1);
  } else {
  }

  /* check for CBF format file */

  *ierr = (cbf_read_file (cbf, in, MSG_DIGEST));
  if(*ierr == 1){
    return(1);
    } 

  /* Get the image identifier */


  cbf_failnez (cbf_rewind_datablock (cbf))

  if (cbf_find_category    (cbf, "diffrn_frame_data") !=0 )

    dimension[0] = *dim_1;
  dimension[1] = *dim_2;


  while (cbf_find_nextrow (cbf, array_id) == 0)
    
      {

	cbf_failnez (cbf_find_column      (cbf, "index"))
	  cbf_failnez (cbf_get_integervalue (cbf, &index))
	  i = index;
	cbf_failnez (cbf_find_column      (cbf, "precedence"))
	  cbf_failnez (cbf_get_integervalue (cbf, &index))
	  if (index >= 1 && index <= 2)
	    {
	      cbf_int[i-1] = index;
	    }
	
	cbf_failnez (cbf_find_column (cbf, "dimension"))
	  cbf_failnez (cbf_get_integervalue (cbf, &dimension[i-1]))
	  cbf_int[i+1] = dimension[i-1];
	
	cbf_failnez (cbf_find_column      (cbf, "direction"))
	  cbf_failnez (cbf_get_value (cbf, &direction[i-1]))
	  
	  strcpy(&cbf_char[i-1][0],direction[i-1]);
	strcat(&cbf_char[i-1][0],"\0");
	
	cbf_failnez (cbf_find_column (cbf, "array_id"))
	  
	  }
    if (dimension [0] == 0 || dimension [1] == 0)
      
      exit (1);
    
    /* Radiation source */
    
    if  (cbf_find_category (cbf, "diffrn_source") == 0) {
      cbf_failnez (cbf_find_column (cbf, "source"))
	cbf_failnez (cbf_get_value   (cbf, &radiation))
	strcpy(&cbf_char[4][0],radiation);
      strcat(&cbf_char[4][0],"\0");
    }
    else {
      strcpy(&cbf_char[4][0],"unspecified radiation source\0");
    }
    i++;
    

  
    /* Detector */

    if (cbf_find_category (cbf, "diffrn_detector") == 0) {
      cbf_failnez (cbf_find_column    (cbf, "id"))
	cbf_failnez(cbf_get_value (cbf, &detector))
	strcpy(&cbf_char[5][0],detector);
      strcat(&cbf_char[5][0],"\0");
    }
    else {
      strcpy(&cbf_char[5][0],"unspecified detector\0");
    }
    i++;
    
    
    /* Pixel size(s) */
    
    if (cbf_find_category    (cbf, "array_element_size") == 0 ){
      cbf_failnez (cbf_find_column    (cbf, "index"))
	cbf_failnez (cbf_get_integervalue (cbf, &index))
	cbf_failnez (cbf_find_column      (cbf, "size"))
	cbf_failnez (cbf_get_doublevalue        (cbf, &pixel_size))
	cbf_double[2] = pixel_size*1000.0;
      
      cbf_failnez(cbf_next_row (cbf))
	cbf_failnez (cbf_find_column    (cbf, "index"))
	cbf_failnez (cbf_get_integervalue (cbf, &index))
	cbf_failnez (cbf_find_column      (cbf, "size"))
	cbf_failnez (cbf_get_doublevalue        (cbf, &pixel_size))
	cbf_double[3] = pixel_size*1000.0;
    }
    else {
      cbf_double[2] = -999.0;
      cbf_double[3] = -999.0;
    }
    /* other pixel information */
    
    if (cbf_find_category  (cbf, "array_intensities") == 0) {
      cbf_failnez (cbf_find_column  (cbf, "gain"))
	cbf_failnez (cbf_get_doublevalue (cbf, &gain))
	cbf_double[4] = gain;
      cbf_failnez (cbf_find_column  (cbf, "overload"))
	cbf_failnez (cbf_get_integervalue (cbf, &overload))
	cbf_int4[0] = overload;
    }
    else {
      cbf_double[4] = -999.0;
      cbf_int4[0] = -9999;
    }
    

  /* Create the new image handle only if *modeop .ne. 2 , otherwise free the 
     file handle (which also closes the file?) */
    if(*modeop == 2){
      cbf_failnez (cbf_free_handle (cbf))
	return 0;
    } 
    else {
      cbf_img = img_make_handle ();
      
      img_set_dimensions (cbf_img, dimension [0], dimension [1]);

  /* Find the binary data */
  
      cbf_failnez (cbf_find_category (cbf, "array_data"))
	cbf_failnez (cbf_find_column   (cbf, "array_id"))
	cbf_failnez (cbf_find_row      (cbf, array_id))
	cbf_failnez (cbf_find_column   (cbf, "data"))
	
	/* Read the binary data */
	//	cbf_failnez (cbf_get_integerarray (cbf, &id, &img_pixel (cbf_img, 0, 0), sizeof(int), 1, img_rows (cbf_img) * img_columns (cbf_img), &nelem_read))
      cbf_failnez (cbf_get_integerarray (cbf, &id, img_pixelptr (cbf_img, 0, 0), sizeof(int), 1, img_rows (cbf_img) * img_columns (cbf_img), &nelem_read))

	/* Free the cbf */
	
	cbf_failnez (cbf_free_handle (cbf))
	
	b = clock ();
      
      /* copy the pixels into an array (short) so that the FORTRAN 
	 side can cope */
      
      /* first we need to decide which way round the image should be read into the 
	 array. This depends on things like the precedence and direction - these 
         are not currently kbnow for Pilatus so we'll just set them to some value */

	dirsta[0] = img_columns(cbf_img)-1;
	dirend[0] = -1;
	dirinc[0] = -1;
	dirsta[1] = 0;
	dirend[1] = img_rows(cbf_img);
	dirinc[1] = 1;
      if (cbf_int[0]==1){
	first = 1;
	second = 0;
      }
      else{
	first = 0;
	second = 1;
      }
      
      colrow = 0;
      //    for (column = dirsta[first]; column != dirend[first]; column = column+dirinc[first])
      //      for (row = dirsta[second]; row != dirend[second];  row=row+dirinc[second])
      for (column = img_columns(cbf_img) - 1; column != 0 ; column--){
	for (row = 0; row < img_rows(cbf_img);  row++){
	  //
	  // superpacking from Andrew
	  //
	  tmp_pixel = img_pixel (cbf_img, column, row);
	  if (tmp_pixel <= 32767) {
	    image[colrow] = tmp_pixel;
	  } 
	  else if (tmp_pixel <131072) {
	    image[colrow] = -((tmp_pixel+8)/16) + 2045;
	    //	    fprintf(stderr,"Pilatus: %d, mosflm: %d\n",tmp_pixel,
	    //    image[colrow]);
	    // fflush(stderr);	  }
	  }
	  else if (tmp_pixel < 3538000) {
	    image[colrow] = -((tmp_pixel+64)/128) - 5124;
	    //	    fprintf(stderr,"pilatus: %d, mosflm: %d\n",tmp_pixel,
	    //	    image[colrow]);
	    //	    fflush(stderr);
	  }
	  else {
	    image[colrow] = -32766;
	  }
	  colrow++;
	}
      }
      
      /*    fprintf (stderr, " Time to read the CBF image: %.3fs\n", 
	  ((b - a) * 1.0) / CLOCKS_PER_SEC);
    */
      
      
      /* Free the images */
      
      img_free_handle (cbf_img);
      
    }
    /* Success */

  return 0;
}
