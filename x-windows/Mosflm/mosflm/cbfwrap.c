/* $Id */
/**********************************************************************
 * cbfwrap -- jiffy to read a CBF file for MOSFLM                      *
 *                                                                    *
 * Version 0.1                                                        *
 *                                                                    *
 *             Harry Powell (harry@mrc-lmb.cam.ac.uk)                 *
 * developed from makecbf by                                          *
 *             Paul Ellis (ellis@ssrl.slac.stanford.edu)              *
 **********************************************************************/
  
/**********************************************************************
 *                                 NOTICE                             *
 * Creative endeavors depend on the lively exchange of ideas. There   *
 * are laws and customs which establish rights and responsibilities   *
 * for authors and the users of what authors create.  This notice     *
 * is not intended to prevent you from using the software and         *
 * documents in this package, but to ensure that there are no         *
 * misunderstandings about terms and conditions of such use.          *
 *                                                                    *
 * Please read the following notice carefully.  If you do not         *
 * understand any portion of this notice, please seek appropriate     *
 * professional legal advice before making use of the software and    *
 * documents included in this software package.  In addition to       *
 * whatever other steps you may be obliged to take to respect the     *
 * intellectual property rights of the various parties involved, if   *
 * you do make use of the software and documents in this package,     *
 * please give credit where credit is due by citing this package,     *
 * its authors and the URL or other source from which you obtained    *
 * it, or equivalent primary references in the literature with the    *
 * same authors.                                                      *
 *                                                                    *
 * Some of the software and documents included within this software   *
 * package are the intellectual property of various parties, and      *
 * placement in this package does not in any way imply that any       *
 * such rights have in any way been waived or diminished.             *
 *                                                                    *
 * With respect to any software or documents for which a copyright    *
 * exists, ALL RIGHTS ARE RESERVED TO THE OWNERS OF SUCH COPYRIGHT.   *
 *                                                                    *
 * Even though the authors of the various documents and software      *
 * found here have made a good faith effort to ensure that the        *
 * documents are correct and that the software performs according     *
 * to its documentation, and we would greatly appreciate hearing of   *
 * any problems you may encounter, the programs and documents any     *
 * files created by the programs are provided **AS IS** without any   *
 * warranty as to correctness, merchantability or fitness for any     *
 * particular or general use.                                         *
 *                                                                    *
 * THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF    *
 * PROGRAMS OR DOCUMENTS OR ANY FILE OR FILES CREATED BY USE OF THE   *
 * PROGRAMS OR DOCUMENTS LIES SOLELY WITH THE USERS OF THE PROGRAMS   *
 * OR DOCUMENTS OR FILE OR FILES AND NOT WITH AUTHORS OF THE          *
 * PROGRAMS OR DOCUMENTS.                                             *
 **********************************************************************/
 
/**********************************************************************
 *                             The IUCr Policy                        *
 *                                    on                              *
 *     the Use of the Crystallographic Information File (CIF)         *
 *                                                                    *
 * The Crystallographic Information File (Hall, Allen & Brown,        *
 * 1991) is, as of January 1992, the recommended method for           *
 * submitting publications to Acta Crystallographica Section C. The   *
 * International Union of Crystallography holds the Copyright on      *
 * the CIF, and has applied for Patents on the STAR File syntax       *
 * which is the basis for the CIF format.                             *
 *                                                                    *
 * It is a principal objective of the IUCr to promote the use of      *
 * CIF for the exchange and storage of scientific data. The IUCr's    *
 * sponsorship of the CIF development was motivated by its            *
 * responsibility to its scientific journals, which set the           *
 * standards in crystallographic publishing. The IUCr intends that    *
 * CIFs will be used increasingly for electronic submission of        *
 * manuscripts to these journals in future. The IUCr recognises       *
 * that, if the CIF and the STAR File are to be adopted as a means    *
 * for universal data exchange, the syntax of these files must be     *
 * strictly and uniformly adhered to. Even small deviations from      *
 * the syntax would ultimately cause the demise of the universal      *
 * file concept. Through its Copyrights and Patents the IUCr has      *
 * taken the steps needed to ensure strict conformance with this      *
 * syntax.                                                            *
 *                                                                    *
 * The IUCr policy on the use of the CIF and STAR File processes is   *
 * as follows:                                                        *
 * _________________________________________________________________  *
 *                                                                    *
 *  * 1 CIFs and STAR Files may be generated, stored or transmitted,  *
 *    without permission or charge, provided their purpose is not     *
 *    specifically for profit or commercial gain, and provided that   *
 *    the published syntax is strictly adhered to.                    *
 *  * 2 Computer software may be developed for use with CIFs or STAR  *
 *    files, without permission or charge, provided it is distributed *
 *    in the public domain. This condition also applies to software   *
 *    for which a charge is made, provided that its primary function  *
 *    is for use with files that satisfy condition 1 and that it is   *
 *    distributed as a minor component of a larger package of         *
 *    software.                                                       *
 *  * 3 Permission will be granted for the use of CIFs and STAR Files *
 *    for specific commercial purposes (such as databases or network  *
 *    exchange processes), and for the distribution of commercial     *
 *    CIF/STAR software, on written application to the IUCr Executive *
 *    Secretary, 2 Abbey Square, Chester CH1 2HU, England. The        *
 *    nature, terms and duration of the licences granted will be      *
 *    determined by the IUCr Executive and Finance Committees.        *
 *                                                                    *
 * _________________________________________________________________  *
 *                                                                    *
 * In summary, the IUCr wishes to promote the use of the STAR File    *
 * concepts as a standard universal data file. It will insist on      *
 * strict compliance with the published syntax for all                *
 * applications. To assist with this compliance, the IUCr provides    *
 * public domain software for checking the logical integrity of a     *
 * CIF, and for validating the data name definitions contained        *
 * within a CIF. Detailed information on this software, and the       *
 * associated dictionaries, may be obtained from the IUCr Office at   *
 * 5 Abbey Square, Chester CH1 2HU, England.                          *
 **********************************************************************/

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
 {int err; \
  err = (x); \
  if (err) { \
    fprintf(stderr,"\nCBFlib fatal error %x \n",err); \
    fprintf(stderr,"caused by call " #x "\n"); \
   } \
 }
/*   local_exit(-1); \ */

int cbfwrap_(int *ierr, int cbf_int[16], long cbf_int4[16], 
	     double cbf_double[16],char cbf_char[16][24],char *argv,
	     short image[64000000],int *modeop)
{
  FILE *in, *out;

  clock_t a,b;

  img_handle img, cbf_img;

  cbf_handle cbf;
  cbf_detector this_detector;
  cbf_goniometer this_goniometer;

  int id, index, iindex;
  int tmp_pixel;

  int i,ii,jjj,colrow, first,second,dirsta[2],dirend[2],dirinc[2];

  unsigned int column, row;

  size_t nelem_read;

  double pixel_size, gain, wavelength, distance, phi_start, phi_range;
  double beam_centre[4];
  double mostimea,mostimeb;
  int overload, dimension[2], precedence ;

  const char *detector,*oscillation_axis,*axis_id;

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
  // ctimer_(&mostimea);
                                           
  /* Create the cbf handle for the image file*/

  cbf_failnez (cbf_make_handle (&cbf))


    /* Read the file */

  in = fopen (argv, "rb");
  if (!in)
  {
    fprintf (stderr, " Couldn't open the CBF file %s\n", argv);
    *ierr = 1;
    return (1);
  }

  /* check for CBF format file */

  *ierr = (cbf_read_file (cbf, in, MSG_DIGEST));
  if(*ierr == 1){
    return(1);
    } 

  /* Get the image identifier */


  cbf_failnez (cbf_rewind_datablock (cbf))

  if (cbf_find_category    (cbf, "diffrn_frame_data") !=0 )
  cbf_failnez (cbf_find_category    (cbf, "diffrn_data_frame"))

  cbf_failnez (cbf_find_column      (cbf, "array_id"))
    
  cbf_failnez (cbf_get_value        (cbf, &array_id))

  /* Get the image dimensions (second dimension = fast, first = slow) */


  cbf_failnez (cbf_find_category    (cbf, "array_structure_list"))

  cbf_failnez (cbf_rewind_row       (cbf))

  cbf_failnez (cbf_find_column      (cbf, "array_id"))

  dimension[0] = dimension[1] = 0;

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

  
  /* Oscillation axis */

  if (cbf_find_category (cbf, "diffrn_measurement_axis") == 0) {
    /*
     *	We can find the scan axis id in the diffrn_measurement_axis category, and look
     *	for that axis in the diffrn_scan_axis category.
     */
    cbf_failnez (cbf_find_column    (cbf, "axis_id"))
      cbf_failnez (cbf_get_value      (cbf, &axis_id))
      
      if (cbf_find_category (cbf, "diffrn_scan_axis") == 0)
	{
	  cbf_failnez (cbf_find_column    (cbf, "axis_id"))
	    do {
	      cbf_failnez (cbf_get_value (cbf, &oscillation_axis))
		if(0 != strcmp(axis_id, oscillation_axis))
		  continue;
	      cbf_failnez (cbf_find_column    (cbf, "angle_start"))
		cbf_failnez (cbf_get_doublevalue (cbf, &phi_start))
		cbf_failnez (cbf_find_column    (cbf, "angle_range"))
		cbf_failnez (cbf_get_doublevalue (cbf, &phi_range))
		cbf_double[5] = phi_start;
	      cbf_double[6] = phi_start + phi_range;
	      cbf_double[7] = phi_range;
	      break;
	    } while(0 == cbf_next_row(cbf));
	}
    
    //    cbf_failnez (cbf_find_column    (cbf, "axis_id"))
    //    cbf_failnez(cbf_get_value (cbf, &oscillation_axis))
    //      strcpy(&cbf_char[3][0],oscillation_axis);
    //    strcat(&cbf_char[3][0],"\0");
  }
  else {	
    /*
     *	No diffrn_measurement_axis category; we will have to try to guess looking at diffrn_scan_axis
     */
    
    if (cbf_find_category (cbf, "diffrn_scan_axis") == 0)
      {
	cbf_failnez (cbf_find_column    (cbf, "axis_id"))
	  do {
	    cbf_failnez (cbf_get_value (cbf, &oscillation_axis))
	      /* if the axis_id has phi or PHI in its name, then we're assuming it's the rotation axis */
	      if(NULL == strstr("phi", oscillation_axis) && NULL == strstr("PHI", oscillation_axis))
		continue;
	    cbf_failnez (cbf_find_column    (cbf, "angle_start"))
	      cbf_failnez (cbf_get_doublevalue (cbf, &phi_start))
	      cbf_failnez (cbf_find_column    (cbf, "angle_range"))
	      cbf_failnez (cbf_get_doublevalue (cbf, &phi_range))
	      cbf_double[5] = phi_start;
	    cbf_double[6] = phi_start + phi_range;
	    cbf_double[7] = phi_range;
	    break;
	  } while(0 == cbf_next_row(cbf));
      }
  }
  
  if(-999.0 == cbf_double[5])
    {
      strcpy(&cbf_char[3][0],"oscillation axis unspecified.\n");
      fprintf(stderr,"oscillation axis unspecified.\n");
    }
  else {
    
    strcpy(&cbf_char[3][0],oscillation_axis);
    strcat(&cbf_char[3][0],"\0");
  }

i++;

  /* hrp 18.01.2007 - new simple way to read header info, using CBFlib 0.7.6.1 + img.c/img.h*/

  cbf_get_wavelength (cbf, &wavelength);
  cbf_double[0] = wavelength;
  cbf_construct_detector (cbf, &this_detector, 0);

  /* b_c[0],[1] slow and fast changing direction in pixels, [2],[3] in mm. */
  
  cbf_get_beam_center (this_detector, &beam_centre [0], &beam_centre [1], &beam_centre [2], &beam_centre [3]);
  
  cbf_double[9] = beam_centre [2];
  cbf_double[10] = -beam_centre [3];
  cbf_get_detector_distance (this_detector, &distance);
  cbf_double[1] = distance;
     
  /* scan axis, angle and size */

  cbf_get_axis_setting (cbf, 0,oscillation_axis, &phi_start, &phi_range);
  cbf_double[5] = phi_start;
  cbf_double[7] = phi_range;
  cbf_double[6] = phi_range + phi_start;
  
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

  /* Polarization of the incident radiation - daft definition at the moment
     so this code is serving only as a place-holder

     if (cbf_find_category  (cbf, "diffrn_radiation_polarizn_ratio") == 0) {
     cbf_failnez (cbf_get_doublevalue (cbf, &polarrat))
     cbf_double[8] = polarrat;
     cbf_failnez (cbf_find_column  (cbf, "polarization_collimation"))
     cbf_failnez (cbf_get_value (cbf, &polarcoll))
     strcpy(&cbf_char[6][0],polarcoll);
     }
     else {
     cbf_double[8] = -999.0;
     strcpy(cbf_char[6][0],"unspecified collimation");
     }
  */


  /* Create the new image handle only if *modeop .ne. 2 , otherwise free the 
     file handle (which also closes the file? do we need the fclose(in)?) */
  /* having an fclose here breaks the code! */
  if(*modeop == 2){
    cbf_failnez (cbf_free_handle (cbf))
      /*      fclose(in); */
    return 0;
  } 
  else {
    cbf_img = img_make_handle ();
  
    if(strcmp(detector,"Pilatus6M") == 0){
      img_set_dimensions (cbf_img, dimension [1], dimension [0]);
    } else {
      img_set_dimensions (cbf_img, dimension [0], dimension [1]);
    }
  

  /* Find the binary data */
  
    cbf_failnez (cbf_find_category (cbf, "array_data"))
      cbf_failnez (cbf_find_column   (cbf, "array_id"))
      cbf_failnez (cbf_find_row      (cbf, array_id))
      cbf_failnez (cbf_find_column   (cbf, "data"))
    
  /* Read the binary data */

      cbf_failnez (cbf_get_integerarray (cbf, &id, img_pixelptr (cbf_img, 0, 0), sizeof(int), 1, img_rows (cbf_img) * img_columns (cbf_img), &nelem_read))

                   /* Free the cbf */
      
      cbf_failnez (cbf_free_handle (cbf))
    
      b = clock ();
  
  /* copy the pixels into an array (short) so that the FORTRAN 
     side can cope */

  /* first we need to decide which way round the image should be read into the 
     array. This depends on things like the precedence and direction */
    if(strcmp (direction[0],"increasing") == 0){
      dirsta[0] = img_columns(cbf_img)-1;
      dirend[0] = -1;
      dirinc[0] = -1;
    }
    else {
      dirsta[0] = 0;
      dirend[0] = img_columns(cbf_img);
      dirinc[0] = 1;
    }
    if(strcmp (direction[1],"increasing") == 0){
      dirsta[1] = 0;
      dirend[1] = img_rows(cbf_img);
      dirinc[1] = 1;
    }
    else{
      dirsta[1] = img_rows(cbf_img)-1;
      dirend[1] = -1;
      dirinc[1] = -1;
    }
    if (cbf_int[0]==1){
      first = 1;
      second = 0;
    }
    else{
      first = 0;
      second = 1;
    }

    colrow = 0;
    
    if(strcmp(detector,"Pilatus6M") == 0){

	 for (column = img_columns(cbf_img) - 1; column != 0 ; column--){
	   for (row = 0; row < img_rows(cbf_img);  row++){
	     /*
	      superpacking from Andrew
	     */
	     tmp_pixel = img_pixel (cbf_img, column, row);
	     if (tmp_pixel <= 32767) {
	       image[colrow] = tmp_pixel;
	     } 
	     else if (tmp_pixel <131072) {
	       image[colrow] = -((tmp_pixel+8)/16) + 2045;
	     }
	     else if (tmp_pixel < 3538000) {
	       image[colrow] = -((tmp_pixel+64)/128) - 5124;
	     }
	     else {
	       image[colrow] = -32766;
	     }
	     colrow++;
	   }
	 }
       } else {
	 for (column = img_columns(cbf_img) - 1; column != 0 ; column--){
	   for (row = 0; row < img_rows(cbf_img);  row++){
	     
	     image[colrow] = img_pixel (cbf_img, column, row);
	     tmp_pixel = image[colrow];
	     if(tmp_pixel < 0){
	       /* conversion from openods.f, around line 1050 */
	       image[colrow] = -(tmp_pixel + 65540)/8;
	     }
	     colrow++;
	   }
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

int local_exit(int status) {
  /* exit(status); */
  return (status);    /* to avoid warning messages */
}
