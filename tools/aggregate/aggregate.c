#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include "mpi.h"
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>

void check_error(int rank, int error_code, int line_number)
{
   char error_string[MPI_MAX_ERROR_STRING];
   int length_of_error_string;

   MPI_Error_string(error_code, error_string, &length_of_error_string);
   fprintf(stderr, "%3d: line %d: %s\n", rank, line_number, error_string);
}

int main(int argc, char **argv)
{
	int ret, Myid,nranks;
	int bufsize;

    MPI_Init( &argc, &argv );
	ret = MPI_Comm_rank(MPI_COMM_WORLD, &Myid);
	if(ret) { check_error(Myid,ret,__LINE__); }
	ret = MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	if(ret) { check_error(Myid,ret,__LINE__); }

	char *dir = argv[1];
	char *filename_pattern = argv[2]; // "alpout.%08d.unw"
	char *output_filename = argv[3];
	char filename[80],f2[80];
	sprintf(f2,filename_pattern,Myid);
	sprintf(filename,"%s/%s",dir,f2);

#ifdef DEBUG
	fprintf(stderr,"%d: using input file %s\n", Myid, filename);
#endif
	int fp = open(filename, O_RDONLY, S_IRUSR | S_IWUSR);


	// stat file to get size
	struct stat s;
	fstat(fp,&s);
	bufsize=s.st_size;

	// read file
	char *buf = (char *)malloc(sizeof(char)*bufsize);
	read(fp,buf,bufsize);
	//fprintf(stderr,"%d: data:%d\n", Myid, buf);
	close(fp);

	// scatter bufsizes
	long long *bufsizes = (long long *)malloc(sizeof(long long)*nranks);
	ret = MPI_Allgather(&s.st_size, 1, MPI_LONG_LONG,
                  bufsizes, 1, MPI_LONG_LONG, 
                  MPI_COMM_WORLD);
	if(ret) { check_error(Myid,ret,__LINE__); }
 
	// compute prefix sum to get offset into file
	long long offset=0; 
	for(int i=1;i<=Myid;i++)
	{
		offset += bufsizes[i-1];
	}

	// write data to rank-specific location
	MPI_File fh;
	MPI_Status status;
	ret = MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if(ret) { check_error(Myid,ret,__LINE__); }
	ret = MPI_File_set_size(fh,0);
	if(ret) { check_error(Myid,ret,__LINE__); }
#ifdef DEBUG
	fprintf(stderr,"%d: Should write %d chars at offset %d\n", Myid, bufsize, offset);
#endif
	ret = MPI_File_write_at_all(fh, offset, buf, bufsize, MPI_CHAR, &status);
	if(ret) { check_error(Myid,ret,__LINE__); }


    int count;
    ret = MPI_Get_elements(&status, MPI_CHAR, &count);
    if(ret != MPI_SUCCESS)
      fprintf(stderr,"MPI_Get_elements error %d\n", ret);
    if(count != bufsize*sizeof(char))
    {
      fprintf(stderr, "Did not write the same number of bytes as requested\n");
      abort();
    }
#ifdef DEBUG
    else
      fprintf(stdout, "Wrote %d bytes\n", count);
#endif


	MPI_File_sync( fh ) ; 
	MPI_Barrier( MPI_COMM_WORLD ) ; 
	ret = MPI_File_close(&fh);
	if(ret) { check_error(Myid,ret,__LINE__); }

	MPI_Finalize();
}
