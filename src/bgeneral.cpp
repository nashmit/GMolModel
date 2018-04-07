#include "bgeneral.hpp"

/**************************************
 * 		General Functions             *
 **************************************/


double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}
 


//
/*
 *  Aminoacids 3 letter notation to 1 letter notation
 */
char aa321 ( const char *aa )// converteste aa din cod de 3 litere in cod de o litera
{
if 	( (!(strcmp(aa, "ALA")) ) || (!(strcmp(aa, "ala")) ) )
	return 'A';
else if ( (!(strcmp(aa, "ARG")) ) || (!(strcmp(aa, "arg")) ) )
	return 'R';
else if ( (!(strcmp(aa, "ASN")) ) || (!(strcmp(aa, "asn")) ) )
	return 'N';
else if ( (!(strcmp(aa, "ASP")) ) || (!(strcmp(aa, "asp")) ) )
	return 'D';
else if ( (!(strcmp(aa, "CYS")) ) || (!(strcmp(aa, "cys")) ) )
	return 'C';
else if ( (!(strcmp(aa, "GLN")) ) || (!(strcmp(aa, "gln")) ) )
	return 'Q';
else if ( (!(strcmp(aa, "GLU")) ) || (!(strcmp(aa, "glu")) ) )
	return 'E';
else if ( (!(strcmp(aa, "GLY")) ) || (!(strcmp(aa, "gly")) ) )
	return 'G';
else if ( (!(strcmp(aa, "HIS")) ) || (!(strcmp(aa, "his")) ) )
	return 'H';
else if ( (!(strcmp(aa, "ILE")) ) || (!(strcmp(aa, "ile")) ) )
	return 'I';
else if ( (!(strcmp(aa, "LEU")) ) || (!(strcmp(aa, "leu")) ) )
	return 'L';
else if ( (!(strcmp(aa, "LYS")) ) || (!(strcmp(aa, "lys")) ) )
	return 'K';
else if ( (!(strcmp(aa, "MET")) ) || (!(strcmp(aa, "met")) ) )
	return 'M';
else if ( (!(strcmp(aa, "PHE")) ) || (!(strcmp(aa, "phe")) ) )
	return 'F';
else if ( (!(strcmp(aa, "PRO")) ) || (!(strcmp(aa, "pro")) ) )
	return 'P';
else if ( (!(strcmp(aa, "SER")) ) || (!(strcmp(aa, "ser")) ) )
	return 'S';
else if ( (!(strcmp(aa, "THR")) ) || (!(strcmp(aa, "thr")) ) )
	return 'T';
else if ( (!(strcmp(aa, "TRP")) ) || (!(strcmp(aa, "trp")) ) )
	return 'W';
else if ( (!(strcmp(aa, "TYR")) ) || (!(strcmp(aa, "tyr")) ) )
	return 'Y';
else if ( (!(strcmp(aa, "VAL")) ) || (!(strcmp(aa, "val")) ) )
	return 'V';
else return 'x';
}

/*
 * aa321 inverse
 */
char aa123 (char *dest, char aa)
{
	if(!dest) return 1;
	if((aa == 'A') || (aa == 'a')){
		strcpy(dest, "ALA");
	}
	else if((aa == 'R') || (aa == 'r')){
		strcpy(dest, "ARG");
	}
	else if((aa == 'N') || (aa == 'n')){
		strcpy(dest, "ASN");
	}
	else if((aa == 'D') || (aa == 'd')){
		strcpy(dest, "ASP");
	}
	else if((aa == 'C') || (aa == 'c')){
		strcpy(dest, "CYS");
	}
	else if((aa == 'Q') || (aa == 'q')){
		strcpy(dest, "GLN");
	}
	else if((aa == 'E') || (aa == 'e')){
		strcpy(dest, "GLU");
	}
	else if((aa == 'G') || (aa == 'g')){
		strcpy(dest, "GLY");
	}
	else if((aa == 'H') || (aa == 'h')){
		strcpy(dest, "HIS");
	}
	else if((aa == 'I') || (aa == 'i')){
		strcpy(dest, "ILE");
	}
	else if((aa == 'L') || (aa == 'l')){
		strcpy(dest, "LEU");
	}
	else if((aa == 'K') || (aa == 'k')){
		strcpy(dest, "LYS");
	}
	else if((aa == 'M') || (aa == 'm')){
		strcpy(dest, "MET");
	}
	else if((aa == 'F') || (aa == 'f')){
		strcpy(dest, "PHE");
	}
	else if((aa == 'P') || (aa == 'p')){
		strcpy(dest, "PRO");
	}
	else if((aa == 'S') || (aa == 's')){
		strcpy(dest, "SER");
	}
	else if((aa == 'T') || (aa == 't')){
		strcpy(dest, "THR");
	}
	else if((aa == 'W') || (aa == 'w')){
		strcpy(dest, "TRP");
	}
	else if((aa == 'Y') || (aa == 'y')){
		strcpy(dest, "TYR");
	}
	else if((aa == 'V') || (aa == 'v')){
		strcpy(dest, "VAL");
	}
	else strcpy(dest, "xxx");

	return 0;
}

/*
 * aa321 inverse
 */
char aa123 (std::string& dest, char aa)
{
	if((aa == 'A') || (aa == 'a')){
		dest = "ALA";
	}
	else if((aa == 'R') || (aa == 'r')){
		dest = "ARG";
	}
	else if((aa == 'N') || (aa == 'n')){
		dest = "ASN";
	}
	else if((aa == 'D') || (aa == 'd')){
		dest = "ASP";
	}
	else if((aa == 'C') || (aa == 'c')){
		dest = "CYS";
	}
	else if((aa == 'Q') || (aa == 'q')){
		dest = "GLN";
	}
	else if((aa == 'E') || (aa == 'e')){
		dest = "GLU";
	}
	else if((aa == 'G') || (aa == 'g')){
		dest = "GLY";
	}
	else if((aa == 'H') || (aa == 'h')){
		dest = "HIS";
	}
	else if((aa == 'I') || (aa == 'i')){
		dest = "ILE";
	}
	else if((aa == 'L') || (aa == 'l')){
		dest = "LEU";
	}
	else if((aa == 'K') || (aa == 'k')){
		dest = "LYS";
	}
	else if((aa == 'M') || (aa == 'm')){
		dest = "MET";
	}
	else if((aa == 'F') || (aa == 'f')){
		dest = "PHE";
	}
	else if((aa == 'P') || (aa == 'p')){
		dest = "PRO";
	}
	else if((aa == 'S') || (aa == 's')){
		dest = "SER";
	}
	else if((aa == 'T') || (aa == 't')){
		dest = "THR";
	}
	else if((aa == 'W') || (aa == 'w')){
		dest = "TRP";
	}
	else if((aa == 'Y') || (aa == 'y')){
		dest = "TYR";
	}
	else if((aa == 'V') || (aa == 'v')){
		dest = "VAL";
	}
	else dest = "xxx";

	return 0;
}


/*
 *  If not, puts '/' at the end of RESULT
 */
void bCheck_path ( char *result, const char *path )
{
	int i;
	strcpy ( result, path );
	i = strlen ( result );
	if ( result[i-1] != '/' )	result[i] = '/';
	result[i+1] = '\0';	
}

/*
 * Puts START-STOP SRC substring into DEST
 * Character is counted from 0
 * Pointers have to be allocated
 */
//indexarea e de la 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//inclusiv caracterul cu indexul stop => '<='
//ai grija: pointerii trebuie sa fie alocati!!!!!!!!!
int bExtractStr ( char *dest, const char *src, int start, int stop )
{
        int i,j;
        for ( i=0,j=start;j<=stop;i++,j++ )     dest[i] = src[j];
        dest[i] = '\0';
        return 0;
}

/*
 * Tolower bExtractStr
 */
int bExtractCaseStr ( char *dest, const char *src, int start, int stop )
{
        int i,j;
        for ( i=0,j=start;j<=stop;i++,j++ )     dest[i] = tolower(src[j]);
        dest[i] = '\0';
        return 0;
}

/*
 * Modifies all characters in a string to NULL
 * Returns how many characters has modified
 */
int bZeroStr(char *dest)
{
	if (dest){
		int len = strlen(dest);
		int i; 
		for(i=0; i<=len; i++){
			dest[i] = '\0';
		}
		return i;
	}
	return 0;
}

/*
 * Modifies NO_ELEM elements in an array to NULL
 * Returns how many elements has modified
 */
int bZeroCharArray(char *array, int no_elem)
{
	if(!array){
		printf("bZeroArray(): NULL pointer passed as argument\n");
		exit(1);	
	}
	int i;
	for(i=0; i<no_elem; i++){
		array[i] = '\0';
	}
	return i;
}

/*
 * Awk Substr
 */
int bSubstr (char *dest, const char *src, int start, int no_chars)
{
	int src_len=strlen(src);
	if((start >= 0) && (no_chars >= 0)){
		int stop=0, k=0, t=0;;
		stop = start + no_chars;
		for(k=start;k<stop;k++,t++){
			if(k>src_len) break;
			dest[t] =  src[k];
		}
		return 0;
	}
	return 1;
}


/*
 * SEQRES line to 1-letter code aa
 * ?????????????????????
 * ?????????????????????
 *
int seqres321(char *dest, const char *src){
	int i;
	char *buffer = new char[4];
	char ter[1];
	for(i=19;i<=70;i+=4){
		bSubstr(buffer,src,i,3);
		ter[1]=aa321(buffer);
		strcat(dest, ter);
	}
	delete buffer;
	return 0;
}
*/

/*
 * Decimal prefix of zeros to limit
 */
string decimal_prefix(double inp_no, long int limit)
{
  if(abs(limit) > 1E+20){
    return string("");
  }
  if((inp_no > (limit/10)) || (inp_no < 0)){
    return string("");
  }

  ostringstream ss;

  int digits = 0;
  long int num = limit;
  while(num){
    num /= 10;
    ++digits;
  }

  string prefix(digits, '\0');
  if(int(inp_no) == 0){
      ss << (limit/10);
      prefix = ss.str();
  }
  for(long int varlim = 10; varlim <= limit; varlim *= 10){
    if((int(inp_no) >= int(varlim/10)) && (inp_no < varlim)){
      ss << (limit/varlim);
      prefix = ss.str();
    } 
  }

  return prefix.substr(1, prefix.length());
}

bool AreSame(double a, double b, double EPSILON)
{
    return fabs(a - b) < EPSILON;
}

/*
 * Convert spatial maatrix (Mat< 2, 2, Mat33 >) to 6x6 matrix (Mat<6,6>)
 * Replaces inf and nan with zeros
 */
bool SpatialMat2Mat66(SimTK::SpatialMat SM, SimTK::Mat<6,6>& M66)
{
    int ti = -1; // 6x6
    int tj = -1; // 6x6
    for(int i = 0; i < 2; i++){ // i for SpatialMatrix
            for(int k = 0; k < 3; k++){ // k for 3x3
        for(int j = 0; j < 2; j++){ // j for SpatialMatrix
                for(int l = 0; l < 3; l++){ // l for 3x3
                    tj++;
                    tj = tj % 6;
                    if(!tj){ti++;}
                    if(std::isinf(SM[i][j][k][l]) || std::isnan(SM[i][j][k][l])){
                        M66[ti][tj] = 0.0;
                    }else{
                        M66[ti][tj] = SM[i][j][k][l];
                    }
                }
            }
        }
    }
    return true;
}

/*
 * Compute numerical matrix inverse with Eigen if possible
 */
bool NumericalInverse(SimTK::Matrix M, SimTK::Matrix& MInv, int nrows, int ncols)
{
    assert((nrows == ncols) && "Robo::NumericalInverse: need a square matrix");
    Eigen::MatrixXd EiM(nrows, ncols);
    Eigen::MatrixXd EiMInv(nrows, ncols);

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            EiM(i, j) = M[i][j];
        }
    }
    
    EiMInv = EiM.inverse();

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            MInv[i][j] = EiMInv(i, j);
        }
    }
    return true;
    
}

/*
 * Compute numerical left inverse with Eigen if possible
 */
bool NumericalLeftInverse(SimTK::Matrix M, SimTK::Matrix& MLeftInv, int nrows, int ncols)
{
    Eigen::MatrixXd EiM(nrows, ncols);
    Eigen::MatrixXd EiMLeftInv(ncols, nrows);

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            EiM(i, j) = M[i][j];
        }
    }
    
    EiMLeftInv = (EiM.transpose() * EiM).inverse() * EiM.transpose();

    for(int i = 0; i < ncols; i++){
        for(int j = 0; j < nrows; j++){
            MLeftInv[i][j] = EiMLeftInv(i, j);
        }
    }
    return true;
    
}

/*
 * Compute numerical right inverse with Eigen if possible
 */
bool NumericalRightInverse(SimTK::Matrix M, SimTK::Matrix& MRightInv, int nrows, int ncols)
{
    Eigen::MatrixXd EiM(nrows, ncols);
    Eigen::MatrixXd EiMRightInv(ncols, nrows);

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            EiM(i, j) = M[i][j];
        }
    }
    
    EiMRightInv = EiM.transpose() * (EiM * EiM.transpose()).inverse();

    for(int i = 0; i < ncols; i++){
        for(int j = 0; j < nrows; j++){
            MRightInv[i][j] = EiMRightInv(i, j);
        }
    }
    return true;
    
}

/*
 * Convert spatial vector to 6-dim vector
 */
SimTK::Vector& SOA_SpatialVec2Vector(SimTK::SpatialVec in, SimTK::Vector& out)
{
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            out[i*3 + j] = in[i][j];
        }
    }
    return out;
}

/*
 * Get the block corresponding to a body from an H-like matrix.
 * Body index "which" starts from 0, 0 being the Ground body.
 */
SimTK::Matrix& SOA_GetHstarLikeElement(SimTK::Matrix inMatrix, int which, SimTK::Matrix& outMatrix)
{
    if(which == 0){ // Ground
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i][j];
            }
        }
    }else if(which == 1){ // First body
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i + 6][j];
            }
        }
    }else{
        for(int i = 0; i < 6; i++){
            outMatrix[i][0] = inMatrix[(which * 6) + i][4 + which];
        }
    }
    return outMatrix;
}

/*
 * Get the block corresponding to a body from an H-like matrix.
 * Body index "which" starts from 0, 0 being the Ground body.
 */
SimTK::Matrix& SOA_GetHLikeElement(SimTK::Matrix inMatrix, int which, SimTK::Matrix& outMatrix)
{
    if(which == 0){ // Ground
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i][j];
            }
        }
    }else if(which == 1){ // First body
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i][j + 6];
            }
        }
    }else{
        for(int i = 0; i < 6; i++){
            outMatrix[0][i] = inMatrix[4 + which][(which * 6) + i];
        }
    }
    return outMatrix;
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Matrix M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat33 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat44 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat66 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Spatial Matrix
 */
void PrintSpatialMat(SimTK::SpatialMat M, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    std::cout << M[i][j][k][l] << " ";
                }
                std::cout << "  ";
            }
        }
        std::cout << std::endl;
    }
}

/*
 * Print Spatial Vector
 */
void PrintSpatialVec(SimTK::SpatialVec M, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places) << std::fixed;
    for(int i = 0; i < 2; i++){
        for(int k = 0; k < 3; k++){
            std::cout << M[i][k] << " ";
        }
    }
    std::cout << std::endl;
}

/*
 * Print Big Vector separated by spaces
 */
void PrintBigMat(SimTK::Vector V, int nrows, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places) << std::fixed;
    for(int i = 0; i < nrows; i++){
            std::cout << V[i] << " ";
    }
    std::cout << std::endl;
}



