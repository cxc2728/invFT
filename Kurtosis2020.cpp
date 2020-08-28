#define _CRT_SECURE_NO_WARNINGS

//  Project Title: Kurtosis Calculator
#include < iostream >
#include < fstream >
#include < string >
#include < io.h >
#include < dos.h >
#include < conio.h >
#include < stdlib.h >
#include < sstream >
#include < stdio.h >
#include < iomanip >
#include < istream >
#include < math.h >

using namespace std;

#define SCALE 255; 

class Kurtosis2D2020 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data { 

		double **Signal; // pointer to the matrix entry

		double **Kurtosis;

		double **average;

		double **m4;

		double **m2;

	}*pointer; // pointer to the matrices

public:

	Kurtosis2D2020(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	void save();

	~Kurtosis2D2020() { } // destructor

};

void Kurtosis2D2020::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;

	 pointer->Signal = new double*[this->n2];

	 pointer->Kurtosis = new double*[this->n2];

	 pointer->average = new double*[this->n2];

	 pointer->m4 = new double*[this->n2];

	 pointer->m2 = new double*[this->n2];


	 for( int v=0; v < this->n2; v++ ) { // (1)

		 pointer->Signal[v] = new double[this->n1];

		 pointer->Kurtosis[v] = new double[this->n1];

		 pointer->average[v] = new double[this->n1];

		 pointer->m4[v] = new double[this->n1];

		 pointer->m2[v] = new double[this->n1];


	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)

			pointer->Signal[v][f] = (double)0.0;

			pointer->Kurtosis[v][f] = (double)0.0;

			pointer->average[v][f] = (double)0.0;

			pointer->m4[v][f] = (double)0.0;

			pointer->m2[v][f] = (double)0.0;

			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


void Kurtosis2D2020::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	

	sprintf(outputFile, "%s","Kurtosis.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->Kurtosis[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","average.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->average[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","m4.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->m4[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","m2.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->m2[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	
} // saveImages



int main ( int argc, char * argv[] ) {

	char outputFile[128]="Kurtosis2D2020.log";

	FILE * savedata;

if (argc < 5) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;			
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the size of the square mask (integer not greater than 20)" << endl;
				 exit(0); }

else { // run the program (begin)

	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	char imageFileName[128];

	int maskSize = atoi(argv[4]);


	if ( maskSize > 20 ) { // if test maskSize
	
		std::cout << "Please make sure that the square maskSize is smaller than 20 pixels" << endl;
	
		fprintf(savedata,"%s\n", "Please make sure that the square maskSize is smaller than 20 pixels");

		fclose(savedata);
	
		exit(0);
	
	} // if test maskSize

	sprintf(imageFileName, "%s", argv[1]);

	std::cout << endl;
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The square mask size is: " << atoi(argv[4]) << endl;


	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The square mask size is:  ", maskSize);

	Kurtosis2D2020 Kurtosis(n1,n2);

	Kurtosis.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		Kurtosis.~Kurtosis2D2020();
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
		
		// at each iteration of the two for loops
		// the program reads the pixel value from the
		// file containing the image and 
		fread(&number,sizeof(double),1,pf);

		Kurtosis.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl; 

	/// Calculate Kurtosis of the image (begins)
	// scale Signal (begin) 
	double max=-Kurtosis.pointer->Signal[0][0];
	double min=Kurtosis.pointer->Signal[0][0];

	for (int i1=0; i1<n2; i1++)
	{ 
			for (int i2=0; i2<n1; i2++)
			{

				if( Kurtosis.pointer->Signal[i1][i2] > (double)max ) 
			
					max = (double)Kurtosis.pointer->Signal[i1][i2];
              
				if( Kurtosis.pointer->Signal[i1][i2] < (double)min ) 
			
					min = (double)Kurtosis.pointer->Signal[i1][i2];
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1<n2; i1++)
	{ 
			for (int i2=0; i2<n1; i2++)
			{

				if ( max == min ) Kurtosis.pointer->Signal[i1][i2] = (double)0.0;

				else { 
					
				Kurtosis.pointer->Signal[i1][i2] = (double) fabs( (min - (double)Kurtosis.pointer->Signal[i1][i2]) / (min - max) );

				Kurtosis.pointer->Signal[i1][i2] *= SCALE;

				// the absolute vale is needed otherwise the minimum is -0 which is nonsense in computers. The value -0 was discovered
				// printing to the screen. Mathematically, the scaling formula gives values in between [0, SCALE] and not [-0, SCALE].
				// Therefore the nonsense is removed taking the absolute value.

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)
	
	
	int n7 = ((int)maskSize);  
	int n8 = ((int)maskSize);

	// pad image memory allocation (begins)
	struct pad_data {

		double **pad_Image; // pointer to the matrix entry 

	}*pad_pointer; // pointer to the matrices

	pad_pointer = new pad_data;

	pad_pointer->pad_Image = new double*[n2+n8];

	for( int v=0; v < n2+n8; v++ ) { // (1)
		 
	pad_pointer->pad_Image[v] = new double[n1+n7];

	  } // (1) allocate struct 'pad_data' (end)

	// (2) initialize (begin)
	for( int v=0; v < n2+n8; v++ ) { // (a)

	    for( int f=0; f < n1+n7 ; f++ ) { // (b)
		 
			pad_pointer->pad_Image[v][f] = (double)0.0;

			 } //(b)

		 } //(a)

	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)
		 
			pad_pointer->pad_Image[v+(int)maskSize/2][f+(int)maskSize/2] = (double)Kurtosis.pointer->Signal[v][f];

			 } //(b)

		 } //(a)
	// (2) initialize (end)

	// pad image memory allocation (ends)

	/// calculate average (begins)
	for (int i = 0; i < (int)n2; i++) { // n1

        for (int j = 0; j < (int)n1; j++) { // n2

			double sum = (double)0.0;

	for (int pp = 0; pp < n8; pp++) { // n7

           for (int qq = 0; qq < n7; qq++) { // n8
			
					sum += (double)pad_pointer->pad_Image[i+pp][j+qq];
          
					} // n8 dim
        
				}  // n7 dim 

			Kurtosis.pointer->average[i][j] = ((double)sum / n7*n8); 

		} // n2 dim
        
	}  // n1 dim 
	/// calculate average (ends)

	
	
	/// calculate Kurtosis (begins)
	for (int i = 0; i < (int)n2; i++) { // n1

        for (int j = 0; j < (int)n1; j++) { // n2

			Kurtosis.pointer->m4[i][j] = (double) 0.0;
			
			Kurtosis.pointer->m2[i][j] = (double) 0.0;

	for (int pp = 0; pp < n8; pp++) { // n7

           for (int qq = 0; qq < n7; qq++) { // n8
			
					Kurtosis.pointer->m4[i][j] += (double) pow( ((double) pad_pointer->pad_Image[i+pp][j+qq] - 
						                                         (double) Kurtosis.pointer->average[i][j]), 4.0);

					Kurtosis.pointer->m2[i][j] += (double) pow( ((double) pad_pointer->pad_Image[i+pp][j+qq] - 
						                                         (double) Kurtosis.pointer->average[i][j]), 2.0);
        
					} // n8 dim
        
				}  // n7 dim 

			if ( (double)Kurtosis.pointer->m2[i][j] != (double) 0.0 )
			
				 Kurtosis.pointer->Kurtosis[i][j] = ( ((double)Kurtosis.pointer->m4[i][j]) / 
				                                      ((double)Kurtosis.pointer->m2[i][j] * Kurtosis.pointer->m2[i][j])); 
			
			else Kurtosis.pointer->Kurtosis[i][j] = (double) 0.0;


		} // n2 dim
        
	}  // n1 dim 
	/// calculate Kurtosis (ends)

	std::cout << "Kurtosis of the image calculated" << endl;
	// calculate Kurtosis (end)

	char PadimageFileName[128];
	sprintf(PadimageFileName, "%s%s", "Pad-", argv[1]);

	// save padded image (begins)
	if ((savedata = fopen(PadimageFileName,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < n2+n8; v++ ) { // (a)

		for( int f=0; f < n1+n7; f++ ) 
	
		fwrite(&pad_pointer->pad_Image[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	// save padded image (ends)


	Kurtosis.save();

	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);

	delete Kurtosis.pointer;
	Kurtosis.~Kurtosis2D2020();
	} // processing (end)

	} // run the program (end)

	
	return 0;
} // end of main 