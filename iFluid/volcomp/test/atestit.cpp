#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "volumecompfnc.h"

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

void adjustSurface(int num_cells, int *gmax_cell, double (*c_corners)[3],int num_points, double (*ps)[3])
{
  int i,j;
  double px,py,pz,dxmin,dxmax,dymin,dymax,dzmin,dzmax;
  double tol=1e-02;

  for (i=0; i < num_points; i++)
  {
    px = ps[i][0];
    py = ps[i][1];
    pz = ps[i][2];
    for (j=0; j< num_cells; j++)
      {
        double cell_bnds[6];
        getCellCorners(j, cell_bnds, gmax_cell, c_corners);      
	
	dxmin = fabs(px-cell_bnds[0])/(cell_bnds[1]-cell_bnds[0]);
        dxmax = fabs(px-cell_bnds[1])/(cell_bnds[1]-cell_bnds[0]);
        dymin = fabs(py-cell_bnds[2])/(cell_bnds[3]-cell_bnds[2]);
        dymax = fabs(py-cell_bnds[3])/(cell_bnds[3]-cell_bnds[2]);
        dzmin = fabs(pz-cell_bnds[4])/(cell_bnds[5]-cell_bnds[4]);
        dzmax = fabs(pz-cell_bnds[5])/(cell_bnds[5]-cell_bnds[4]);

        if (dxmin < tol)
          {
           ps[i][0] = cell_bnds[0];
          }
       else if (dxmax < tol)
          {
           ps[i][0] = cell_bnds[1];
	  }
        if (dymin < tol)
          {
	    ps[i][1] = cell_bnds[2];
	  }
        else if (dymax < tol)
          {
	    ps[i][1] = cell_bnds[3];
	  }
        if (dzmin < tol)
          {
	    ps[i][2] = cell_bnds[4];
	  }
        else if (dzmax < tol)
          {
	    ps[i][2] = cell_bnds[5];
	  }
      }
  }
}


int main(int argc, char **argv) {

  /* INPUT DECLARATION*/

  int i; 
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
 
  FILE *infile_cell = fopen(argv[1], "r"); 
  FILE *infile_intfc = fopen(argv[2], "r");
 

  /* Read Intfc Information */
  fscanf(infile_intfc, "%d", &num_points);
  ps = new double [num_points][3];

  for (i = 0; i < num_points; i++) {
      fscanf(infile_intfc, "%lf", &ps[i][0]);
      fscanf(infile_intfc, "%lf", &ps[i][1]);
      fscanf(infile_intfc, "%lf", &ps[i][2]);
  }


  fscanf(infile_intfc, "%d", &num_tris);
  tris = new int [num_tris][3];

  comp_index = new int [num_tris];
  for (i = 0; i < num_tris; i++) {
      fscanf(infile_intfc, "%d", &tris[i][0]);
      fscanf(infile_intfc, "%d", &tris[i][1]);
      fscanf(infile_intfc, "%d", &tris[i][2]);
      fscanf(infile_intfc, "%d", &comp_index[i] );
  }

  fscanf(infile_intfc, "%d", &num_comps);
  comp_list = new int [num_comps];
  pos_map = new int [num_comps];
  neg_map = new int [num_comps];
  pos_vol = new double [num_comps];
  neg_vol = new double [num_comps];

  for (i = 0; i < num_comps; i++) {
      fscanf(infile_intfc, "%d", &comp_list[i]);
      pos_map[i] = neg_map[i] = comp_list[i];
      pos_vol[i] = neg_vol[i] = 0.0;
  }

  /***************************/

  /* Read Cell Information */

  fscanf(infile_cell, "%d", &num_cells);
  c_centers = new double [num_cells][3];
  face_type = new int [num_cells][6];
  face_sign = new unsigned char [num_cells][6]; 
  
  for (i = 0; i < num_cells; i++)
  {
      fscanf(infile_cell, "%lf", &c_centers[i][0]);
      fscanf(infile_cell, "%lf", &c_centers[i][1]);
      fscanf(infile_cell, "%lf", &c_centers[i][2]);

      for (int j = 0; j < 6; j++) {
	  face_type[i][j] = -2;
	  face_sign[i][j] = 1;
      }
  }

  fscanf(infile_cell, "%lf", &(grid_spacing[0]));
  fscanf(infile_cell, "%lf", &(grid_spacing[1]));
  fscanf(infile_cell, "%lf", &(grid_spacing[2]));

  fscanf(infile_cell, "%d", &num_corners);
    c_corners = new double [num_corners][3];
      for (i = 0; i < num_corners; i++)
      {
	fscanf(infile_cell, "%lf", &c_corners[i][0]);
      fscanf(infile_cell, "%lf", &c_corners[i][1]);
      fscanf(infile_cell, "%lf", &c_corners[i][2]);  
      }
      
  fscanf(infile_cell, "%d", &(gmax_cell[0]));
  fscanf(infile_cell, "%d", &(gmax_cell[1]));
  fscanf(infile_cell, "%d", &(gmax_cell[2]));

      

  fclose(infile_intfc);
  fclose(infile_cell);

  printf ("Finished Reading Input !\n");
  
  printf ("Starting Preprocessing \n");
  adjustSurface(num_cells,gmax_cell,c_corners,num_points,ps);
  printf ("Finished Preprocessing \n");

  /*Call the entry function*/

  for (int k = 0; k < num_cells; k++){
      double cell_bnds[6];
      getCellCorners(k, cell_bnds, gmax_cell, c_corners);
      /* printf("Cell Bounds = %g %g %g %g %g %g\n",cell_bnds[0],cell_bnds[1],cell_bnds[2],cell_bnds[3],cell_bnds[4],cell_bnds[5]);*/
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
