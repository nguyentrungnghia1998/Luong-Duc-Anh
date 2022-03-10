#include "lbm_D2Q9.h"


int D2Q9_output_vtk(int nout, double **rho, double **u, double **v, double **type)
{
	int km = 1;
	int		i,j,k;
	char	filename[128];
	FILE	*fp;
	unsigned int array_size;
	unsigned long int offset=0;
	//short  num16; // Int16 2byte
	float  val32; // Float32 4byte

	sprintf(filename,"./field%06d.vtr",nout);
	fp=fopen(filename,"wb");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
	fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",im,jm,km);
	fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",im,jm,km);
	fprintf(fp,"    <PointData>\n");
	fprintf(fp,"    </PointData>\n");
	fprintf(fp,"    <CellData>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(im)*(jm)*(km);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(im)*(jm)*(km);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(im)*(jm)*(km);
	fprintf(fp,"    </CellData>\n");
	fprintf(fp,"    <Coordinates>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(im+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(jm+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(km+1);
	fprintf(fp,"    </Coordinates>\n");
	fprintf(fp,"  </Piece>\n");
	fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"  <AppendedData encoding=\"raw\">");
	fprintf(fp,"_");

    // Density (cell)
	array_size=1*4*(im)*(jm)*(km);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<km+1;k++){
		for(j=1;j<jm+1;j++){
			for(i=1;i<im+1;i++){
				val32=(float)rho[i][j]; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}
    // Velocity (cell)
	array_size=3*4*(im)*(jm)*(km);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<km+1;k++){
		for(j=1;j<jm+1;j++){
			for(i=1;i<im+1;i++){
				val32=(float)u[i][j]; fwrite(&val32,sizeof(float),1,fp);
				val32=(float)v[i][j]; fwrite(&val32,sizeof(float),1,fp);
				val32=0.0; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}
    // CellType (cell)
	array_size=1*4*(im)*(jm)*(km);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<km+1;k++){
		for(j=1;j<jm+1;j++){
			for(i=1;i<im+1;i++){
				val32=(float)type[i][j]; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}
	// Coordinates (vertices)
	array_size=1*4*(im+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(i=0;i<im+1;i++){ val32=(float)(i*dr); fwrite(&val32,sizeof(float),1,fp); }
	array_size=1*4*(jm+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(j=0;j<jm+1;j++){ val32=(float)(j*dr); fwrite(&val32,sizeof(float),1,fp); }
	array_size=1*4*(km+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<km+1;k++){ val32=(float)(k*dr); fwrite(&val32,sizeof(float),1,fp); }
	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);

	nout++;

	return nout;
}