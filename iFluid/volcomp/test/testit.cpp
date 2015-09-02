#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "volumecompfnc.h"

typedef int L3ARRAY [3];

bool hardware_is_little_endian()
{
    	short int word = 0x0001;
    	char *byte = (char *) &word;
    	if(byte[0] == 1) 
	    return true;
    	else
            return false;
}

int endian_int_swap(int i)
{
        int j;
        int const size = sizeof(int);
        union
	{
	    int i; 
            unsigned char b[16];
        } dat1, dat2;
        
	dat1.i = i;
        for(j = 0; j < size; ++j)
	    dat2.b[j] = dat1.b[size-1-j];

	return dat2.i;
}

double endian_double_swap(double f)
{
        int j;
        const int size = sizeof(double);
        union
        {
            double f;
            unsigned char b[16];
        } dat1, dat2;

        dat1.f = f;
        for(j = 0; j < size; ++j)
   	    dat2.b[j] = dat1.b[size-1-j];

        return dat2.f;
}
void getCellCorners(int index, double *cell_bnds, int *gmax_cell, double (*c_corners)[3])
{
    int icoords[3];
    int id = index;

    for (int d = 0; d < 3; d++)
    {
	int G = gmax_cell[d];
	icoords[d] = id % G;
	id = (id - icoords[d])/G;
    }

    int xmin_index = icoords[0];
    int xmax_index = icoords[0] + 1;
    int ymin_index = icoords[1];
    int ymax_index = icoords[1] + 1;
    int zmin_index = icoords[2];
    int zmax_index = icoords[2] + 1;

    int Lcorner_index = (xmin_index + ymin_index*(gmax_cell[0]+1) + zmin_index*(gmax_cell[0]+1)*(gmax_cell[1]+1));
    int Ucorner_index = (xmax_index + ymax_index*(gmax_cell[0]+1) + zmax_index*(gmax_cell[0]+1)*(gmax_cell[1]+1));

    cell_bnds[0] = c_corners[Lcorner_index][0]; //xmin
    cell_bnds[1] = c_corners[Ucorner_index][0]; //xmax

    cell_bnds[2] = c_corners[Lcorner_index][1]; //ymin
    cell_bnds[3] = c_corners[Ucorner_index][1]; //ymax

    cell_bnds[4] = c_corners[Lcorner_index][2]; //zmin
    cell_bnds[5] = c_corners[Ucorner_index][2]; //zmax
}


int main(int argc, char **argv) {

  /* INPUT DECLARATION*/

  int i;
  int ival[4];
  double val[3];
  int num_points, num_tris, num_comps, num_cells, num_corners;
  double grid_spacing[3];
  double (*c_centers)[3];
  double (*c_corners)[3];
  int gmax_cell[3];
  int    (*face_type)[6];
  unsigned char (*face_sign)[6];

 
  int *comp_list;
  int *pos_map;
  int *neg_map;
  double *pos_vol;
  double *neg_vol;

  double (*ps)[3];
  int    (*tris)[3];
  int *comp_index;
  
 
  FILE *infile_cell = fopen(argv[1], "rb");
  FILE *infile_intfc = fopen(argv[2], "rb");

  if (hardware_is_little_endian())
  {
      fread(ival, sizeof(int), 1, infile_intfc);
      num_points = endian_int_swap(ival[0]);
      ps = new double [num_points][3];

      for (i = 0; i < num_points; i++)
      {
	  fread(val, sizeof(double), 3, infile_intfc);
	  ps[i][0] = endian_double_swap(val[0]);
	  ps[i][1] = endian_double_swap(val[1]);
	  ps[i][2] = endian_double_swap(val[2]);
      }

      fread(ival, sizeof(int), 1, infile_intfc);
      num_tris = endian_int_swap(ival[0]);
      tris = new int [num_tris][3];
      comp_index = new int [num_tris];

      for (i = 0; i < num_tris; i++)
      {
	  fread(ival, sizeof(int), 4, infile_intfc);
	  tris[i][0] = endian_int_swap(ival[0]);
	  tris[i][1] = endian_int_swap(ival[1]);
	  tris[i][2] = endian_int_swap(ival[2]);
	  comp_index[i] = endian_int_swap(ival[3]);
      }

      fread(ival, sizeof(int), 1, infile_intfc);
      num_comps = endian_int_swap(ival[0]);
      comp_list = new int [num_comps];
      pos_map = new int [num_comps];
      neg_map = new int [num_comps];
      pos_vol = new double [num_comps];
      neg_vol = new double [num_comps];

      for (i = 0; i < num_comps; i++)
      {
	  fread(ival, sizeof(int), 1, infile_intfc);
	  comp_list[i] = endian_int_swap(ival[0]);
	  pos_map[i] = neg_map[i] = comp_list[i];
	  pos_vol[i] = neg_vol[i] = 0.0;
      }

      fread(ival, sizeof(int), 1, infile_cell);
      num_cells = endian_int_swap(ival[0]);
      c_centers = new double [num_cells][3];
      face_type = new int [num_cells][6];
      face_sign = new unsigned char [num_cells][6];
      for (i = 0; i < num_cells; i++)
      {
	  fread(val, sizeof(double), 3, infile_cell);
	  c_centers[i][0] = endian_double_swap(val[0]);
	  c_centers[i][1] = endian_double_swap(val[1]);
	  c_centers[i][2] = endian_double_swap(val[2]);
	  for (int j = 0; j < 6; j++) {
	      face_type[i][j] = -2;
	      face_sign[i][j] = 1;
	  }
      }

      fread(val, sizeof(double), 3, infile_cell);
      grid_spacing[0] = endian_double_swap(val[0]);
      grid_spacing[1] = endian_double_swap(val[1]);
      grid_spacing[2] = endian_double_swap(val[2]);

      fread(ival, sizeof(int), 1, infile_cell);
      num_corners = endian_int_swap(ival[0]);
      c_corners = new double [num_corners][3];
      for (i = 0; i < num_corners; i++)
      {
	  fread(val, sizeof(double), 3, infile_cell);
	  c_corners[i][0] = endian_double_swap(val[0]);
	  c_corners[i][1] = endian_double_swap(val[1]);
	  c_corners[i][2] = endian_double_swap(val[2]);
      }

      fread(ival, sizeof(int), 3, infile_cell);
      gmax_cell[0] = endian_int_swap(ival[0]);
      gmax_cell[1] = endian_int_swap(ival[1]);
      gmax_cell[2] = endian_int_swap(ival[2]);
  }	   
  else
  {
      fread(ival, sizeof(int), 1, infile_intfc);
      num_points = ival[0];
      ps = new double [num_points][3];

      for (i = 0; i < num_points; i++)
      {
	  fread(val, sizeof(double), 3, infile_intfc);
	  ps[i][0] = val[0];
	  ps[i][1] = val[1];
	  ps[i][2] = val[2];
      }

      fread(ival, sizeof(int), 1, infile_intfc);
      num_tris = ival[0];
      tris = new int [num_tris][3];
      comp_index = new int [num_tris];

      for (i = 0; i < num_tris; i++)
      {
	  fread(ival, sizeof(int), 4, infile_intfc);
	  tris[i][0] = ival[0];
	  tris[i][1] = ival[1];
	  tris[i][2] = ival[2];
	  comp_index[i] = ival[3];
      }
      fread(ival, sizeof(int), 1, infile_intfc);
      num_comps = ival[0];
      comp_list = new int [num_comps];
      pos_map = new int [num_comps];
      neg_map = new int [num_comps];
      pos_vol = new double [num_comps];
      neg_vol = new double [num_comps];
      for (i = 0; i < num_comps; i++)
      {
	  fread(ival, sizeof(int), 1, infile_intfc);
	  comp_list[i] = ival[0];
	  pos_map[i] = neg_map[i] = comp_list[i];
	  pos_vol[i] = neg_vol[i] = 0.0;
      }


      fread(ival, sizeof(int), 1, infile_cell);
      num_cells = ival[0];
      c_centers = new double [num_cells][3];
      face_type = new int [num_cells][6];
      face_sign = new unsigned char [num_cells][6];
      for (i = 0; i < num_cells; i++)
      {
	  fread(val, sizeof(double), 3, infile_cell);
	  c_centers[i][0] = val[0];
	  c_centers[i][1] = val[1];
	  c_centers[i][2] = val[2];
	  for (int j = 0; j < 6; j++) {
	      face_type[i][j] = -2;
	      face_sign[i][j] = 1;
	  }
      }

      fread(val, sizeof(double), 3, infile_cell);
      grid_spacing[0] = val[0];
      grid_spacing[1] = val[1];
      grid_spacing[2] = val[2];

      fread(ival, sizeof(int), 1, infile_cell);
      num_corners = ival[0];
      c_corners = new double [num_corners][3];
      for (i = 0; i < num_corners; i++)
      {
	  fread(val, sizeof(double), 3, infile_cell);
	  c_corners[i][0] = val[0];
	  c_corners[i][1] = val[1];
	  c_corners[i][2] = val[2];
      }
      fread(ival, sizeof(int), 3, infile_cell);
      gmax_cell[0] = ival[0];
      gmax_cell[1] = ival[1];
      gmax_cell[2] = ival[2];
  }	   

  fclose(infile_intfc);
  fclose(infile_cell);

  printf ("Finished Reading !\n");

  /*Call the entry function*/

  for (int k = 0; k < num_cells; k++){
      double cell_bnds[6];
      getCellCorners(k, cell_bnds, gmax_cell, c_corners);
      volumecompfnc(cell_bnds,face_type[k],face_sign[k], ps[0],num_points,
		    tris[0],num_tris,comp_index, num_comps,comp_list, pos_map, neg_map, pos_vol,neg_vol);
      printf("Finished Cell Number = %d\n",k);
      /*printf("Cell centers = %.6g %.6g %.6g\n",c_centers[k][0],c_centers[k][1],c_centers[k][2]);*/
      printf("Face type = %d %d %d %d %d %d \n",face_type[k][0],face_type[k][1],face_type[k][2],face_type[k][3],face_type[k][4],face_type[k][5]);
      printf("Face sign = %d %d %d %d %d %d \n",face_sign[k][0],face_sign[k][1],face_sign[k][2],face_sign[k][3],face_sign[k][4],face_sign[k][5]);

  }
 
  delete [] ps;
  delete [] tris;

  delete [] c_centers[0];
  delete [] face_type[0];
  delete [] face_sign[0];

  delete [] comp_index;
  delete [] comp_list;
  delete [] pos_map;
  delete [] neg_map;
  delete [] pos_vol;
  delete [] neg_vol;

  return(0);
}
