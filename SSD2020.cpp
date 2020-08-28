#define _CRT_SECURE_NO_WARNINGS

//  Project Title: Sum of Square Differences - SSD Calculator
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

class SSD2D2020 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data { 

		double **ICF; // pointer to the matrix entry

		double **reconICF; // pointer to the matrix entry

		double **SSD; // pointer to the matrix entry



	}*pointer; // pointer to the matrices

public:

	SSD2D2020(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	void save();

	~SSD2D2020() { } // destructor

};

void SSD2D2020::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;

	 pointer->ICF = new double*[this->n2];

	 pointer->reconICF = new double*[this->n2];

	 pointer->SSD = new double*[this->n2];


	 for( int v=0; v < this->n2; v++ ) { // (1)

		 pointer->ICF[v] = new double[this->n1];

		 pointer->reconICF[v] = new double[this->n1];

		 pointer->SSD[v] = new double[this->n1];

	
	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)

			pointer->ICF[v][f] = (double)0.0;

			pointer->reconICF[v][f] = (double)0.0;

			pointer->SSD[v][f] = (double)0.0;


			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


void SSD2D2020::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	

	sprintf(outputFile, "%s","SSD.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->SSD[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	
} // saveImages



int main ( int argc, char * argv[] ) {

	char outputFile[128]="SSD2D2020.log";

	FILE * savedata;

if (argc < 6) { std::cout << endl;
				 std::cout << "Please type the ICF image file name" << endl;
				 std::cout << "Please type the reconICF image file name" << endl;
				 std::cout << "Please make sure that the images format is Analyze 'double': 64 bits real" << endl;			
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the size of the square mask (integer not greater than 20)" << endl;
				 exit(0); }

else { // run the program (begin)

	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	int n1 = atoi(argv[3]);
	int n2 = atoi(argv[4]);

	char imageFileName[128];
	char imageFileName_2[128];

	int maskSize = atoi(argv[5]);


	if ( maskSize > 20 ) { // if test maskSize
	
		std::cout << "Please make sure that the square maskSize is smaller than 20 pixels" << endl;
	
		fprintf(savedata,"%s\n", "Please make sure that the square maskSize is smaller than 20 pixels");

		fclose(savedata);
	
		exit(0);
	
	} // if test maskSize

	sprintf(imageFileName, "%s", argv[1]);
	sprintf(imageFileName_2, "%s", argv[2]);

	std::cout << endl;
	std::cout << "The ICF image file name is: " << imageFileName << endl;
	std::cout << "The reconICF image file name is: " << imageFileName_2 << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[3]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[4]) << endl;
	std::cout << "The square mask size is: " << atoi(argv[5]) << endl;


	fprintf(savedata,"%s%s\n", "The ICF image file name is: " , imageFileName);
	fprintf(savedata,"%s%s\n", "The reconICF image file name is: " , imageFileName_2);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The square mask size is:  ", maskSize);

	SSD2D2020 SSD(n1,n2);

	SSD.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		SSD.~SSD2D2020();
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
		
		// at each iteration of the two for loops
		// the program reads the pixel value from the
		// file containing the image and 
		fread(&number,sizeof(double),1,pf);

		SSD.pointer->ICF[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)


	if ((pf = fopen(imageFileName_2,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName_2 << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName_2 );
		SSD.~SSD2D2020();
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
		
		// at each iteration of the two for loops
		// the program reads the pixel value from the
		// file containing the image and 
		fread(&number,sizeof(double),1,pf);

		SSD.pointer->reconICF[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl; 
	
	int n7 = ((int)maskSize);  
	int n8 = ((int)maskSize);

	// pad image memory allocation (begins)
	struct pad_data {

		double **pad_Image;  // pointer to the matrix entry (ICF)

		double **pad_Image2;  // pointer to the matrix entry (reconICF)

	}*pad_pointer; // pointer to the matrices

	pad_pointer = new pad_data;

	pad_pointer->pad_Image = new double*[n2+n8];

	pad_pointer->pad_Image2 = new double*[n2+n8];

	for( int v=0; v < n2+n8; v++ ) { // (1)
		 
	pad_pointer->pad_Image[v] = new double[n1+n7];

	pad_pointer->pad_Image2[v] = new double[n1+n7];

	  } // (1) allocate struct 'pad_data' (end)

	// (2) initialize (begin)
	for( int v=0; v < n2+n8; v++ ) { // (a)

	    for( int f=0; f < n1+n7 ; f++ ) { // (b)
		 
			pad_pointer->pad_Image[v][f] = (double)0.0;

			pad_pointer->pad_Image2[v][f] = (double)0.0;

			 } //(b)

		 } //(a)

	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)
		 
			pad_pointer->pad_Image[v+(int)maskSize/2][f+(int)maskSize/2] = (double)SSD.pointer->ICF[v][f];

			pad_pointer->pad_Image2[v+(int)maskSize/2][f+(int)maskSize/2] = (double)SSD.pointer->reconICF[v][f];

			 } //(b)

		 } //(a)
	// (2) initialize (end)

	// pad image memory allocation (ends)

	/// calculate sum of square differences - SSD (begins)
	for (int i = 0; i < (int)n2; i++) { // n1

        for (int j = 0; j < (int)n1; j++) { // n2

			double sum = (double)0.0;

	for (int pp = 0; pp < n8; pp++) { // n7

           for (int qq = 0; qq < n7; qq++) { // n8
			
					sum += pow( ((double)pad_pointer->pad_Image[i+pp][j+qq] - (double)pad_pointer->pad_Image2[i+pp][j+qq]), 2.0);
          
					} // n8 dim
        
				}  // n7 dim 

			SSD.pointer->SSD[i][j] = (double)sum; 

		} // n2 dim
        
	}  // n1 dim 
	/// calculate sum of square differences - SSD (ends)

	std::cout << "SSD of the images calculated" << endl;
	// calculate SSD (end)

	char PadimageFileName[128];
	char PadimageFileName2[128];

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


	sprintf(PadimageFileName2, "%s%s", "Pad-", argv[2]);

	// save padded image (begins)
	if ((savedata = fopen(PadimageFileName2,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < n2+n8; v++ ) { // (a)

		for( int f=0; f < n1+n7; f++ ) 
	
		fwrite(&pad_pointer->pad_Image2[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	// save padded image (ends)

	SSD.save();

	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);

	delete SSD.pointer;
	SSD.~SSD2D2020();
	} // processing (end)

	} // run the program (end)

	
	return 0;
} // end of main 