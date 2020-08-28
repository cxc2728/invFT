//  This file contains sourcecode distributed as freeware. 
//  The intellectual property of the sourcecode is shown 
//  here to belong to Carlo Ciulla.

// Disclaimer: 

// The website here named www.sourcecodewebsiteCarloCiulla.com [1] does not intend 
// to convey the meaning of profit making for what pertains to the content
// provided. --->>> Instead, when the content is downloaded, the user(s) are
// kindly invited to donate money to charity organizations involved in 
// helping people in need of food and water. <<<---


// The Novel Re-sampling Locations have been sized to be a fraction of 
// the pixel size. The programs presented here confirm both concepts and 
// implications brought to knowledge through the unifying theory [1].

// Reference:

// [1] Carlo Ciulla "Improved Signal and Image Interpolation in Biomedical Applications: 
// The Case of Magnetic Resonance Imaging (MRI)." Medical Information Science 
// Reference - IGI Global Publisher - March 2009; ISBN: 978 - 160566202 - 2.

//  Project Title: Phase Calculator
#define _CRT_SECURE_NO_WARNINGS

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


#define ft_SCALE 255 
// this constant is necessary in order to control the 
// contrast brightness of the Fourier Transformations 
// (direct and inverse)


using namespace std;

void OnInverseFourierTransformMerge(char icfFileName[], int rcxres, int rcyres, char imageFileName[]);
void OnFourierTransform(char imageFilename[], int rcxres, int rcyres);
void OnInverseFourierTransform(char filename[], int rcyres, int rcxres);
void OnCalculatePhase(int rcxres, int rcyres, char imageFileName[]);
void OnInverseFourierTransformKSpace(char icfFileName[], int rcxres, int rcyres, char imageFileName[]);

class SRE2D2013 {

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

		double **ICF_SRE2D; // pointer to the matrix entry

		double **SRE2D_ictBI; // pointer to the matrix entry

		double **SRE2D_ictAI; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	SRE2D2013(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	void save();

	~SRE2D2013() { } // destructor

};

void SRE2D2013::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;
			
	 pointer->Signal = new double*[this->n2];

	 pointer->ICF_SRE2D = new double*[this->n2];

	 pointer->SRE2D_ictBI = new double*[this->n2];

	 pointer->SRE2D_ictAI = new double*[this->n2];


	 for( int v=0; v < this->n2; v++ ) { // (1)
		 
		 pointer->Signal[v] = new double[this->n1];

		 pointer->ICF_SRE2D[v] = new double[this->n1];

		 pointer->SRE2D_ictBI[v] = new double[this->n1];

		 pointer->SRE2D_ictAI[v] = new double[this->n1];


	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for(int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)
		 
			pointer->Signal[v][f] = (double)0.0;

			pointer->ICF_SRE2D[v][f] = (double)0.0;

			pointer->SRE2D_ictBI[v][f] = (double)0.0;

			pointer->SRE2D_ictAI[v][f] = (double)0.0;

		    } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


void SRE2D2013::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	
	sprintf(outputFile, "%s","Signal.img");

	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->Signal[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","ICF_SRE2D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->ICF_SRE2D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","SRE2D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->SRE2D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","SRE2D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		fwrite(&pointer->SRE2D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

} // saveImages



int main ( int argc, char * argv[] ) {

	char outputFile[128]="Phase2D2020.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 9) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				 std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the X direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the Y direction (double)" << endl;
				 std::cout << "Please enter the XY rotation angle (double)" << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)
	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);
	
	double XPixelSize = atof(argv[4]);
	double YPixelSize = atof(argv[5]);

	double x_misplacement_X = atof(argv[6]);
	double y_misplacement_Y = atof(argv[7]);

	double theta = atof(argv[8]);
	
	char imageFileName[128];
	
	sprintf(imageFileName, "%s", argv[1]);
	
	std::cout << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The XY rotation angle is: " << atof(argv[8]) << endl;
	
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);

    double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
    double misplacement_Y = ((double)1.0 - ( -sin( (double)theta ) + cos( (double)theta ) ) + y_misplacement_Y);

      misplacement_X = ((double)misplacement_X/XPixelSize);
      misplacement_Y = ((double)misplacement_Y/YPixelSize);

	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////


	SRE2D2013 SRE(n1,n2);

	SRE.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
			
		fread(&number,sizeof(double),1,pf);
		
		SRE.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;

	// calculate ICF_SRE2D (begin)
	std::cout << "Compute ICF of the SRE2D model function" << endl;
	
	double k = 0.0, s = 0.0;
	for (int i1=0; i1 < n2-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n1-1; i2++) { // y dim

		  
        k = (double) SRE.pointer->Signal[i1][i2] * misplacement_X * misplacement_Y + 
				   ( misplacement_X * misplacement_Y * misplacement_X / 2.0) * 
				   ( SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2] ) + 
				   ( misplacement_X * misplacement_Y * misplacement_Y / 2.0) * 
				   ( SRE.pointer->Signal[i1][i2+1] - SRE.pointer->Signal[i1][i2] ) +
				   ( misplacement_X * misplacement_X * misplacement_Y * misplacement_Y / 4.0) * 
				   ( SRE.pointer->Signal[i1+1][i2+1] + SRE.pointer->Signal[i1][i2] - SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2+1] );


        s = (double) SRE.pointer->Signal[i1][i2] * misplacement_X * misplacement_Y; 

		SRE.pointer->SRE2D_ictBI[i1][i2] = (double)s;

		SRE.pointer->SRE2D_ictAI[i1][i2] = (double)k;
       
		   if ( (double)k != 0.0 ) SRE.pointer->ICF_SRE2D[i1][i2] = (double)s/k;  
		   else SRE.pointer->ICF_SRE2D[i1][i2] = (double) 0.0;

		} // y dim
        
	}  // x dim


	std::cout << "ICF of the SRE2D model function calculated" << endl;
	// calculate ICF_SRE2D (end)

	// scale ICF_SRE2D (begin) 
	double max=-MAX;
	double min=MAX;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
	
		if( SRE.pointer->ICF_SRE2D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_SRE2D[i1][i2];
              
		if( SRE.pointer->ICF_SRE2D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_SRE2D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_SRE2D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_SRE2D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_SRE2D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_SRE2D (end)


	SRE.save(); // save all of the images


	std::cout << "Please Wait...Now running the Fourier transformation on the Signal" << endl;

	OnFourierTransform(imageFileName, n1, n2);

	OnCalculatePhase(n1, n2, imageFileName);

	std::cout << "Please Wait...Now running the Fourier transformation on the Phase" << endl;
	OnFourierTransform("base-phase.img", n1, n2);
	
	char reconFilename[128];

	sprintf(reconFilename, "%s%s", "mri-phase-", imageFileName);
	OnFourierTransform(reconFilename, n1, n2);

	sprintf(reconFilename, "%s%s", "phase-sum-", imageFileName);
	OnFourierTransform(reconFilename, n1, n2);

	std::cout << "Please Wait...Now running the Fourier transformation on the ICF" << endl;
	OnFourierTransform("ICF_SRE2D.img", n1, n2);

	std::cout << "Inverse FT of the Signal and the ICF_SRE2D... with base phase" << endl;
	OnInverseFourierTransformMerge("ICF_SRE2D.img", n1, n2, imageFileName);

	std::cout << "Inverse FT of the Signal and the ICF_SRE2D...w/o base phase" << endl;
	OnInverseFourierTransformKSpace("ICF_SRE2D.img", n1, n2, imageFileName);
	
	
	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);

	delete SRE.pointer;
	SRE.~SRE2D2013();
	} // processing (end)

	} // run the program (end)

	return 0;
} // end of main 

void OnInverseFourierTransformMerge(char icfFileName[], int rcxres, int rcyres, char imageFileName[])
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;
	
	double phase, readData;

	//2010
	double emittingSource = 1.0; 
	double scale = ((double)rcxres*rcyres*emittingSource); 
	//2010

	FILE * logfile;
	char logfilename[128]="FourierTransformMerge.log";

	FILE *image;
	char imageFilename[256];
	FILE * icf_pf;

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * reconSignal = 0;
	double * icfReal = 0;
	double * icfImaginary = 0;

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now FT Processing...");
    fprintf(logfile,"%s\n", "Now FT Processing...");

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}


	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Imaginary Image data: Exit");
	
		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);

		exit(0);

	}

	
	if ((icfReal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
    {

            fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

            // FIFO memory deallocation method
            free(kSpaceR);
            free(kSpaceI);
            free(reconSignal);
            exit(0);

    }

	if ((icfImaginary = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
    {

            fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

            // FIFO memory deallocation method
            free(kSpaceR);
            free(kSpaceI);
            free(reconSignal);
            free(icfReal);
			exit(0);

    }

} // allocate memory

	
	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", imageFileName);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
     free(icfImaginary);
		
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR+index) = (double) number;

			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "K-SpaceI-", imageFileName);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);
		
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI+index) = (double) number;

			}

		}

		fclose(image);

	}// read data and initialize pointers

	char theFile[128];

	sprintf(theFile,  "%s%s", "K-SpaceR-", icfFileName);

    if ((icf_pf = fopen(theFile,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);

	 exit(0);
	
	} 
	
	else { // read data

	    //Read file
        for (i=0; i<rcyres; i++)
        { ///read signal data
            for (j=0; j<rcxres; j++)
            {

                index = ((j*rcyres)+i);

                fread(&readData,sizeof(double),1,icf_pf);
			
                *(icfReal + index) = (double)readData;

            }
		}

        fclose(icf_pf);

	} ///read data

	
	sprintf(theFile, "%s%s", "K-SpaceI-", icfFileName);

    if ((icf_pf = fopen(theFile,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
     free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);

	 exit(0);
	
	} 
	
	else { // read data

	    //Read file
        for (i=0; i<rcyres; i++)
        { ///read signal data
            for (j=0; j<rcxres; j++)
            {

                index = ((j*rcyres)+i);

                fread(&readData,sizeof(double),1,icf_pf);
			
				*(icfImaginary + index) = (double)readData;
     
            }
		}

        fclose(icf_pf);
	} ///read data

	double real = 0.0, imaginary = 0.0;

	///// Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{
		
	    	dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);
		
	  	    k2 = ((int)(dx*NofXpixels)+dy);

			w = ((j*NofYpixels)+i);

			real = (double) 0.0;
			imaginary = (double) 0.0;

			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);
					
					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					double mriPhase = (double)2.0*pi*atan2( (double) *(kSpaceI+t) , (double) *(kSpaceR+t)  ) / (NofXpixels*NofYpixels);
							
					phase += (double)mriPhase;

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					
					double complexR = (double) cos( (double)phase ) + (double) sin( (double)phase ); 

					double complexI = -(double) sin( (double)phase ) + (double) cos( (double)phase ); 

					
					real += (double) *(icfReal+t) * (double) complexR;

					imaginary -= (double) *(icfImaginary+t) * (double) complexI;

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
			}

		}///process k-space data 

			
			*(reconSignal+w) =  (double) sqrt( ((double) real * real)  + ((double) imaginary * imaginary) );

    		*(reconSignal+w) /= (double)scale;
			
		}
	} ///process k-space data


	// scale *(reconSignal) (begin) 
	double max=-*(reconSignal);
	double min=*(reconSignal);

	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

			w = ((j*NofYpixels)+i);

				if( (double)*(reconSignal+w) > (double)max ) 
			
					max = (double)*(reconSignal+w);
              
				if( (double)*(reconSignal+w) < (double)min ) 
			
					min = (double)*(reconSignal+w);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

				w = ((j*NofYpixels)+i);

				if ( max == min ) (double)*(reconSignal+w);

				else (double)*(reconSignal+w)  =  (double)min - ( ( (double)*(reconSignal+w) * (min - max) ) / (double) ft_SCALE);
			
		} // y dim
        
	}  // x dim
	// scale *(reconSignal+w) (end)

	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	/// save the reconstructed signal (begins)
	sprintf(reconFilename, "%s%s", "recon-Signal-W-base-phase", icfFileName);


    fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Reconstructed Signal Saved");

	fclose (pf);
	} // save data
	/// save the reconstructed signal (ends)
		
}

void OnInverseFourierTransformKSpace(char icfFileName[], int rcxres, int rcyres, char imageFileName[])
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;
	
	double readData;

	//2010
	double emittingSource = 1.0; 
	double scale = ((double)rcxres*rcyres*emittingSource); 
	//2010

	FILE * logfile;
	char logfilename[128]="FourierTransformKSpace.log";

	FILE *image;
	char imageFilename[256];
	FILE * icf_pf;

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * reconSignal = 0;
	double * icfReal = 0;
	double * icfImaginary = 0;

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now FT Processing...");
    fprintf(logfile,"%s\n", "Now FT Processing...");

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}


	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Imaginary Image data: Exit");
	
		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);

		exit(0);

	}

	
	if ((icfReal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
    {

            fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

            // FIFO memory deallocation method
            free(kSpaceR);
            free(kSpaceI);
            free(reconSignal);
            exit(0);

    }

	if ((icfImaginary = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
    {

            fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

            // FIFO memory deallocation method
            free(kSpaceR);
            free(kSpaceI);
            free(reconSignal);
            free(icfReal);
			exit(0);

    }

	} // allocate memory

	
	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", imageFileName);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
     free(icfImaginary);
		
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR+index) = (double) number;
		
			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "K-SpaceI-", imageFileName);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);
		
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI+index) = (double) number;

			}

		}

		fclose(image);

	}// read data and initialize pointers

	char theFile[128];

	sprintf(theFile,  "%s%s", "K-SpaceR-", icfFileName);

    if ((icf_pf = fopen(theFile,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);

	 exit(0);
	
	} 
	
	else { // read data

	    //Read file
        for (i=0; i<rcyres; i++)
        { ///read signal data
            for (j=0; j<rcxres; j++)
            {

                index = ((j*rcyres)+i);

                fread(&readData,sizeof(double),1,icf_pf);
			
                *(icfReal + index) = (double)readData;

            }
		}

        fclose(icf_pf);

	} ///read data

	
	sprintf(theFile, "%s%s", "K-SpaceI-", icfFileName);

    if ((icf_pf = fopen(theFile,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
     free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);

	 exit(0);
	
	} 
	
	else { // read data

	    //Read file
        for (i=0; i<rcyres; i++)
        { ///read signal data
            for (j=0; j<rcxres; j++)
            {

                index = ((j*rcyres)+i);

                fread(&readData,sizeof(double),1,icf_pf);
			
				*(icfImaginary + index) = (double)readData;
     
            }
		}

        fclose(icf_pf);
	} ///read data

	double real = 0.0, imaginary = 0.0;

	///// Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{
		
	    	dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);
		
	  	    k2 = ((int)(dx*NofXpixels)+dy);

			w = ((j*NofYpixels)+i);

			real = (double) 0.0;
			imaginary = (double) 0.0;

			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);

					double mriPhase = (double)2.0*pi*atan2( (double) *(kSpaceI+t) , (double) *(kSpaceR+t)  ) / (NofXpixels*NofYpixels);

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					
					double complexR = (double) cos( (double)mriPhase ) + (double) sin( (double)mriPhase ); 

					double complexI = -(double) sin( (double)mriPhase ) + (double) cos( (double)mriPhase ); 

					
					real += (double) *(icfReal+t) * (double) complexR;

					imaginary -= (double) *(icfImaginary+t) * (double) complexI;

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
			}

		}///process k-space data 

			
			*(reconSignal+w) =  (double) sqrt( ((double) real * real)  + ((double) imaginary * imaginary) );

    		*(reconSignal+w) /= (double)scale;
			
		}
	} ///process k-space data


	// scale *(reconSignal) (begin) 
	double max=-*(reconSignal);
	double min=*(reconSignal);

	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

			w = ((j*NofYpixels)+i);

				if( (double)*(reconSignal+w) > (double)max ) 
			
					max = (double)*(reconSignal+w);
              
				if( (double)*(reconSignal+w) < (double)min ) 
			
					min = (double)*(reconSignal+w);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

				w = ((j*NofYpixels)+i);

				if ( max == min ) (double)*(reconSignal+w);

				else (double)*(reconSignal+w)  =  (double)min - ( ( (double)*(reconSignal+w) * (min - max) ) / (double) ft_SCALE);
			
		} // y dim
        
	}  // x dim
	// scale *(reconSignal+w) (end)

	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	/// save the reconstructed signal (begins)
	sprintf(reconFilename, "%s%s", "recon-Signal-WO-base-phase", icfFileName);


    fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	 free(icfReal);
	 free(icfImaginary);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Reconstructed Signal Saved");

	fclose (pf);
	} // save data
	/// save the reconstructed signal (ends)
		
}

void OnCalculatePhase(int rcxres, int rcyres, char imageFileName[])
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k3, k2, w, t;
	
	double pi = 3.141592;

	//2010
	double emittingSource = 1.0; 
	double scale = ((double)rcxres*rcyres*emittingSource); 
	//2010

	FILE * logfile;
	char logfilename[128]="CalculatePhase.log";

	FILE *image;
	char imageFilename[256];

	double * kSpaceR = 0;
	double * kSpaceI = 0;

	class phase2020 {

		int nx;
	    int ny;

		public:

				struct data {

				double **basephase;
				double **mriphase;
				double **sumphase;

				}*pointer;

		public:

			phase2020(int x, int y) : nx(x), ny(y) { }; // constructor 
	
			void allocatePhase() {  // allocate

					pointer = new data;

					pointer->basephase = new double*[this->ny];
					
					pointer->mriphase = new double*[this->ny];
					
					pointer->sumphase = new double*[this->ny];


					 for( int v=0; v < this->ny; v++ ) { // (1)
			
						pointer->basephase[v] = new double[this->nx];

						pointer->mriphase[v] = new double[this->nx];

						pointer->sumphase[v] = new double[this->nx];

					} // (1) allocate struct 'data' (end)

					 
					 for(int v=0; v < this->ny; v++ )

						for( int f=0; f < this->nx ; f++ ) { // (b)
			
						pointer->basephase[v][f] = (double)0.0;

						pointer->mriphase[v][f] = (double)0.0;

						pointer->sumphase[v][f] = (double)0.0;
			} //(b)

		} // allocate

			~phase2020() { } // destructor
	};

		phase2020 phase(rcxres, rcyres);

		phase.allocatePhase();


	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now Calculating the Phase...");
    fprintf(logfile,"%s\n", "Now Calculating the Phase...");

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}

	} // allocate memory

	
	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", imageFileName);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);

	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "K-SpaceI-", imageFileName);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);

	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI+index) = (double) number;

			}

		}

		fclose(image);

	}// read data and initialize pointers


			///// Phase Calculator (begins) //////
			for (i=0; i<NofYpixels; i++)
			{ ///process k-space data

				for (j=0; j<NofXpixels; j++)
				{
		
	    			dx = ((int) i - NofYpixels/2);
					dy = ((int) j - NofXpixels/2);
		
	  				k2 = ((int)(dx*NofXpixels)+dy);

					w = ((j*NofYpixels)+i);


			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);
					
					double bphase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					phase.pointer->basephase[s][p] = (double) bphase / ((double) NofXpixels*NofYpixels*NofXpixels*NofYpixels);

					phase.pointer->mriphase[s][p] = (double)2.0*pi*atan2( (double) *(kSpaceI+t) , (double) *(kSpaceR+t)  );
		
				}
			}

		
		}
	} ///process k-space data
	
			// sum the MRI phase to the ICF phase (begins)
			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 
					
					phase.pointer->sumphase[s][p] = ( (double) phase.pointer->basephase[s][p] + 
						                              (double) phase.pointer->mriphase[s][p] ); 	
				}
			} 
			// sum the MRI phase to the ICF phase (ends)
			///// Phase Calculator (ends) //////

	char reconFilename[128];
	FILE * pf;
	double savedata;
	

	// save the phase (begins)
	sprintf(reconFilename, "%s", "base-phase.img");


    fprintf(logfile, "%s\t%s\n", "Now Saving Phase Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			savedata = phase.pointer->basephase[i][j]; 
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Phase Signal Saved");

	fclose (pf);
	} // save data
	// save the icf phase (ends)
	
	// save the mri phase (begins)
	sprintf(reconFilename, "%s%s", "mri-phase-", imageFileName);


    fprintf(logfile, "%s\t%s\n", "Now Saving Phase Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			savedata = (double) phase.pointer->mriphase[i][j]; 
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Phase Signal Saved");

	fclose (pf);
	} // save data
	// save the mri phase (ends)

	// save the sum phase (begins)
	sprintf(reconFilename, "%s%s", "phase-sum-", imageFileName);


    fprintf(logfile, "%s\t%s\n", "Now Saving Phase Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			savedata = (double) phase.pointer->sumphase[i][j];
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Phase Signal Saved");

	fclose (pf);
	} // save data
	// save the sum phase (ends)

	delete phase.pointer;
	phase.~phase2020();
}


void OnFourierTransform(char imageFilename[], int rcxres, int rcyres)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;

	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * Signal = 0;

	FILE * logfile;
	
	char logfilename[128]="Fourier-T.log";

  	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

		printf("%s\n %s\n" , "Unable to open log File", "Now Exit");

		exit(0);
	
	} else { // allocate memory 


	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}

	if ((Signal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		exit(0);

	}

	} // allocate memory 

	//// read image data and initialize pointers
	double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);
				
				*(kSpaceR+index) = (double) 0.0;

				*(kSpaceI+index) = (double) 0.0;

			}

		}

	FILE * pf;
	char SignalFilename[128];
	double readData;
	
	sprintf(SignalFilename, "%s", imageFilename);

	if ((pf = fopen(SignalFilename,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // read data


	for (i=0; i<rcyres; i++)
	{ ///read signal data
		for (j=0; j<rcxres; j++)
		{

			index = ((j*rcyres)+i);
          
            fread(&readData,sizeof(double),1,pf);

			*(Signal+index) = (double)readData;

		}
	} ///read signal data

	fprintf(logfile,"%s\n", "Signal Read in DOUBLE (64bits) format");

	fclose (pf);
	} // save data

	double phase, complexR, complexI;
	
	///// Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///calculate k-space data

		for (j=0; j<NofXpixels; j++)
		{

	
			dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);

			k2 = ((int)(dy*NofYpixels)+dx); 

			w = ((j*NofYpixels)+i);

			for (int s=0; s<NofYpixels; s++)
			{ ///calculate k-space data 
				for (int p=0; p<NofXpixels; p++)
				{ 
					

		     		ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);
 
				    k3 = ((int)(ds*NofXpixels)+dp); 

					t = ((p*NofYpixels)+s);

					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					complexR = (double) cos( (double)phase ) + (double) sin( (double)phase ); 

					complexI = -(double) sin( (double)phase ) + (double) cos( (double)phase ); 
					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
				
					*(kSpaceR+w) += (double) *(Signal+t) * (double) complexR;

					*(kSpaceI+w) -= (double) *(Signal+t) * (double) complexI;

			}

		}///calculate k-space data 

			    
		}
	} ///calculate k-space data

	///// Fourier Transform //////
	double savedata = 0.0;
	char FTfilename[128];

	sprintf(FTfilename, "%s%s", "K-SpaceR-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Real) in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");


	 // FIFO memory deallocation method
 	 free(kSpaceR);
 	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(kSpaceR+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Real) Saved");

	fclose (pf);
	} // save data



	sprintf(FTfilename, "%s%s", "K-SpaceI-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Imaginary) in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(kSpaceI+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Imaginary) Saved");

	fclose (pf);
	
	} // save data

	sprintf(FTfilename, "%s%s", "K-SpaceM-", imageFilename);

    fprintf_s(logfile, "%s\t%s\n", "Now Saving K-Space Magnitude of the Signal in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save K-Space Magnitude of the Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data	

		// save a zero image (begin)
		for (int s=0; s<NofYpixels; s++)
		{ 
			for (int p=0; p<NofXpixels; p++)
			{ 

			savedata = (double)0.0;
          
            fwrite(&savedata,sizeof(double),1,pf);

			}
		} // save a zero image (end)

	fclose(pf);
	
	}
		
	if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save K-Space Magnitude of the Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data
		
		// K-Space Magnitude (begin)
		for (int s=0; s<(int)NofYpixels; s++)
		{ 
			for (int p=0; p<(int)NofXpixels; p++)
			{ 
			
		
			index = ((p*NofYpixels)+s);

			savedata = (double) sqrt( (double)*(kSpaceR+index)*(double)*(kSpaceR+index) + 
		   		                      (double)*(kSpaceI+index)*(double)*(kSpaceI+index) );
          
            fwrite(&savedata,sizeof(double),1,pf);
			
		}
	} // K-Space Magnitude (end)

	fprintf_s(logfile,"%s\n", "K-Space Magnitude of the Signal Saved");

	fclose (pf);
	} // save data

	printf("%s\n", "FT Processing Completed");
    fprintf_s(logfile,"%s\n", "FT Processing Completed");

	fclose(logfile);
	
	// FIFO memory deallocation method
	free(kSpaceR);
	free(kSpaceI);
	free(Signal);

}

void OnInverseFourierTransform(char filename[], int rcxres, int rcyres)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;
	
	double phase;

	//2010
	double emittingSource = 0.98; 
	double scale = ((double)rcxres*rcyres*emittingSource); 
	//2010

	FILE * logfile;
	char logfilename[128]="INV-FourierT.log";

	FILE *image;
	char imageFilename[256];

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * reconSignal = 0;

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now INV FT Processing...");
    fprintf(logfile,"%s\n", "Now INV FT Processing...");

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}


	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Imaginary Image data: Exit");
	
		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);

		exit(0);

	}

	} // allocate memory

	
	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", filename);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "K-SpaceI-", filename);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI+index) = (double) number;

			}

		}

		fclose(image);


		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				*(reconSignal+index) = (double)0.0;
					
			}

		}


	}// read data and initialize pointers

	double real = 0.0, imaginary = 0.0;
	
	///// INV Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{
		
	    	dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);
		
	  	    k2 = ((int)(dx*NofXpixels)+dy);

			w = ((j*NofYpixels)+i);

			real = (double)0.0;
			imaginary = (double)0.0;

			
			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);
					
					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					real += ((double) *(kSpaceR+t) * cos( (double) phase)) + ((double) *(kSpaceI+t) * (double) sin((double)phase));

					imaginary += -((double) *(kSpaceR+t) * sin((double)phase)) + ((double) *(kSpaceI+t) * cos((double)phase)); 
					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
			}

		}///process k-space data 

			*(reconSignal+w) =  (double) sqrt( ((double) real * real)  + ((double) imaginary * imaginary) );

			*(reconSignal+w) /= (double)scale;
		}
	} ///process k-space data



	// scale *(reconSignal) (begin) 
	double max=-*(reconSignal);
	double min=*(reconSignal);

	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

				w = ((j*NofYpixels)+i);

				if( *(reconSignal+w) > (double)max ) 
			
					max = (double)*(reconSignal+w);
              
				if( *(reconSignal+w) < (double)min ) 
			
					min = (double)*(reconSignal+w);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

				w = ((j*NofYpixels)+i);

				if ( max == min ) *(reconSignal+w) = (double) 0.0;

				else *(reconSignal+w)  =  (double)min - ( ( (double)*(reconSignal+w) * (min - max) ) / (double) ft_SCALE);
			
		} // y dim
        
	}  // x dim
	// scale *(reconSignal+w) (end)

	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	sprintf(reconFilename, "%s%s", "reconSignal-", filename);


    fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Reconstructed Signal Saved");

	fclose (pf);
	} // save data


    printf("%s\n", "Inverse FT Processing Completed");
    fprintf(logfile,"%s\n", "Inverse FT Processing Completed");

	fclose(logfile);
		
	// FIFO memory deallocation method
	free(kSpaceR);
	free(kSpaceI);
	free(reconSignal);

}