/*
	Francois Belanger 94 245 437
	Genevieve Dostie
*/

#include "MImage.h"
#include "stdio.h"
#include "stdlib.h"
#include <cmath>
#include "gc/GCoptimization.h"
#include <limits>
#include <algorithm>

MImage::MImage(int xs,int ys,int zs)
{
	MXS = 0;
	MYS = 0;
	MZS = 0;
	MImgBuf=NULL;

	if(xs>0 && ys>0 && zs>0)
		MAllocMemory(xs,ys,zs);

	for(int y=0;y<MYS;y++)
		for(int x=0;x<MXS;x++)
			MSetColor(0,x,y);
}
MImage::MImage(void)
{
	MXS = 0;
	MYS = 0;
	MZS = 0;
	MImgBuf = NULL;
}

MImage::MImage(const MImage &copy)
{
	MXS = 0;
	MYS = 0;
	MZS = 0;
	MImgBuf = NULL;
	if (copy.MIsEmpty()) {
		return;
	}

	MAllocMemory(copy.MXS, copy.MYS, copy.MZS);

	for (int y = 0; y<MYS; y++)
		for (int x = 0; x<MXS; x++)
			MImgBuf[x][y] = copy.MImgBuf[x][y];
}

MImage::~MImage(void)
{
	MFreeMemory();
}

/* =================================================================================
====================================================================================
======================                    I/O                 ======================
====================================================================================
====================================================================================*/

/*
Allocates a xs x ys x zs image buffer
*/
void MImage::MAllocMemory(int xs, int ys, int zs)
{
	MFreeMemory();
	if (xs <= 0 || ys <= 0 || zs <= 0) {
		printf("Error !! MImage::MAllocMemory\n");
		return;
	}

	MXS = xs;
	MYS = ys;
	MZS = zs;

	MImgBuf = new RGBPixel*[xs];

	for (int i = 0; i<xs; i++)
		MImgBuf[i] = new RGBPixel[ys];
}

/*
Free the image buffer
*/
void MImage::MFreeMemory(void)
{
	if (MImgBuf == NULL || MXS <= 0 || MYS <= 0 || MZS <= 0) {
		MXS = MYS = MZS = 0;
	}
	else {
		for (int i = 0; i<MXS; i++)
			delete[]MImgBuf[i];
		delete[]MImgBuf;
		MImgBuf = NULL;
		MXS = MYS = MZS = 0;
	}
}

/*
load a pgm/ppm image
*/
bool MImage::MLoadImage(const char *fileName)
{

	if (fileName == NULL) {
		printf("Error!!! MImage::MLoadImage\n");
		return false;
	}
	FILE_FORMAT ff;
	char tmpBuf[500];
	std::ifstream inFile;
	int maxVal, val;
	char valRaw;
	unsigned char color;

	inFile.open(fileName, std::ios::in);
	if (!inFile.is_open()) {
		printf("Error!!! cant open file %s\n", fileName);
		return false;
	}

	inFile.getline(tmpBuf, 500);
	switch (tmpBuf[1]) {
	case '2':
		ff = PGM_ASCII;
		MZS = 1;
		break;
	}

	int nbComm = 0;
	inFile.getline(tmpBuf, 500);
	while (tmpBuf[0] == '#') {
		nbComm++;
		inFile.getline(tmpBuf, 500);
	}
	inFile.close();

	if (ff == PGM_ASCII)
		inFile.open(fileName, std::ios::in);
	else
		inFile.open(fileName, std::ios::in | std::ios::binary);

	inFile.getline(tmpBuf, 500);
	while (nbComm>0) {
		nbComm--;
		inFile.getline(tmpBuf, 500);
	}

	inFile >> MXS;
	inFile >> MYS;
	inFile >> maxVal;

	MAllocMemory(MXS, MYS, MZS);
	switch (ff) {

	case PGM_ASCII:
		for (int y = 0; y<MYS; y++)
			for (int x = 0; x<MXS; x++) {
				inFile >> val;
				MImgBuf[x][y].r = (float)val*255.0 / maxVal;
			}
		break;

	default:
		printf("Error!!! cant read that format. Please use PGM_ASCII\n");
		inFile.close();
		return false;
	}

	inFile.close();
	printf("File %s opened successfully\n", fileName);
	return true;
}

/*
save a pgm/ppm image
*/
bool MImage::MSaveImage(const char *fileName, FILE_FORMAT ff)
{
	if (!fileName) {
		printf("Error!! MImage::MSaveImage\n");
		return false;
	}
	unsigned char val;
	std::ofstream outFile;

	outFile.open(fileName, std::ios::out | std::ios::binary);
	if (!outFile.is_open()) {
		printf("Error!!! cant open file %s\n", fileName);
		return false;
	}
	switch (ff) {
	case PGM_ASCII:
		outFile << "P2" << "\n" << MXS << " " << MYS << "\n" << "255" << "\n";
		for (int y = 0; y<MYS; y++) {
			for (int x = 0; x<MXS; x++) {
				outFile << (unsigned short)(MImgBuf[x][y].r) << " ";
			}
			outFile << "\n";
		}
		break;

	default:
		printf("Error!!! can't write that format. Please use PGM_ASCII\n");
		outFile.close();
		return false;
	}
	outFile.close();
	printf("File %s saved successfully\n", fileName);
	return true;
}



/* =================================================================================
====================================================================================
======================           Point Operations             ======================
====================================================================================
====================================================================================*/


/*
Every pixel with an intensity > 'tvalue' are set to 255.  The other ones are set to 0.
*/
void MImage::MThreshold(float tvalue)
{
	for (int y = 0; y<MYS; y++)
		for (int x = 0; x<MXS; x++) {

			if (MZS>1) {

				if ((MImgBuf[x][y].r + MImgBuf[x][y].g + MImgBuf[x][y].b) / 3.0>tvalue) {
					MImgBuf[x][y].r = 255;
					MImgBuf[x][y].g = 255;
					MImgBuf[x][y].b = 255;
				}
				else {
					MImgBuf[x][y].r = 0;
					MImgBuf[x][y].g = 0;
					MImgBuf[x][y].b = 0;
				}

			}
			else {
				if (MImgBuf[x][y].r>tvalue)
					MImgBuf[x][y].r = 255;
				else
					MImgBuf[x][y].r = 0;
			}
		}
}

/*
Rescale the image between 0 and 255
*/
void MImage::MRescale(void)
{
	float maxR, minR;
	float maxG, minG;
	float maxB, minB;

	maxR = maxG = maxB = -99999;
	minR = minG = minB = 99999;
	int X, Y;
	for (int y = 0; y<MYS; y++)
		for (int x = 0; x<MXS; x++) {

			if (MImgBuf[x][y].r>maxR)
				maxR = MImgBuf[x][y].r;
			if (MImgBuf[x][y].r<minR)
				minR = MImgBuf[x][y].r;

			if (MZS>1) {
				if (MImgBuf[x][y].g>maxG)
					maxG = MImgBuf[x][y].g;
				if (MImgBuf[x][y].b>maxB)
					maxB = MImgBuf[x][y].b;

				if (MImgBuf[x][y].g<minG)
					minG = MImgBuf[x][y].g;
				if (MImgBuf[x][y].b<minB)
					minB = MImgBuf[x][y].b;
			}
		}

	for (int y = 0; y<MYS; y++)
		for (int x = 0; x<MXS; x++) {
			MImgBuf[x][y].r = (MImgBuf[x][y].r - minR)*255.0 / (maxR - minR);

			if (MZS>1) {
				MImgBuf[x][y].g = (MImgBuf[x][y].g - minG)*255.0 / (maxG - minG);
				MImgBuf[x][y].b = (MImgBuf[x][y].b - minB)*255.0 / (maxB - minB);
			}
		}
}




/*
Mean shift filtering

the implementation is inspired of the following paper
D. Comanicu, P. Meer: "Mean shift: A robust approach toward feature space analysis".
IEEE Trans. Pattern Anal. Machine Intell., May 2002.
The resulting filtered image is copied in the current image (this->MImgBuf)
*/
void MImage::MMeanShift(float SpatialBandWidth, float RangeBandWidth, float tolerance)
{

}
/* =================================================================================
====================================================================================
======================            Feature extraction          ======================
====================================================================================
====================================================================================*/

void innondation(MImage &X_s, std::vector<std::vector<int>> &Y_s, int xSeed, int ySeed, float tolerance, float color_target) {
	std::vector<std::pair<int, int>> x_y;
	x_y.push_back(std::pair<int, int>(xSeed + 1, ySeed));
	x_y.push_back(std::pair<int, int>(xSeed - 1, ySeed));
	x_y.push_back(std::pair<int, int>(xSeed, ySeed + 1));
	x_y.push_back(std::pair<int, int>(xSeed, ySeed - 1));

	for (auto i = x_y.begin(); i != x_y.end(); i++) {
		if (i->first >= 0 && i->first < X_s.MXSize() && i->second >= 0 && i->second < X_s.MYSize()) {
			if (Y_s[i->first][i->second] < 0 && std::abs(X_s.MGetColor(i->first, i->second) - color_target) <= tolerance) {
				Y_s[i->first][i->second] = 1;
				innondation(X_s, Y_s, i->first, i->second, tolerance, color_target);
			}
		}
	}
}

/*
Segmentation with magic wand algorithm
(xSeed, ySeed) is where the region starts growing and
tolerance is the criteria to stop the region from growing.
*/
void MImage::MMagicWand(int xSeed, int ySeed, float tolerance)
{
	std::vector<std::vector<int>> map_black_white;
	map_black_white = std::vector<std::vector<int> >(MXS, std::vector<int>(MYS, -1));
	map_black_white[xSeed][ySeed] = 1;
	innondation(*this, map_black_white, xSeed, ySeed, tolerance, MGetColor(xSeed, ySeed));

	for (int y = 0; y < MYSize(); y++)
	{
		for (int x = 0; x < MXSize(); x++)
		{
			if (map_black_white[x][y] > 0) {
				MImgBuf[x][y].r = 255;
			}
			else {
				MImgBuf[x][y].r = 0;
			}
		}
	}
}



/*
N-class segmentation.  Each class is a Gaussian function defined by
a mean, stddev and prior value.
The resulting label Field is copied in the current image (this->MImgBuf)
*/
void MImage::MOptimalThresholding(float *means, float *stddev, float *apriori, int nbClasses)
{
	if (nbClasses <= 1)
	{
		std::cout << "Error : number of classes should be at least 2!" << std::endl;
		return;
	}

	int classRange = floor(255.0f / nbClasses);
	float sqTwoPi = sqrtf(2.0f * M_PI);

	for (int y = 0; y < MYSize(); y++)
	{
		for (int x = 0; x < MXSize(); x++)
		{
			int c = 0;
			float maxProb = 0.0f;
			for (int k = 0; k < nbClasses; k++)
			{
				float diff = powf(MGetColor(x, y) - means[k], 2.0f);
				float prob = apriori[k] * exp2f(-diff / (2.0f * powf(stddev[k], 2.0f))) / (sqTwoPi * stddev[k]);
				if (prob > maxProb)
				{
					maxProb = prob;
					c = k;
				}
			}
			MSetColor(c * classRange, x, y);
		}
	}
}

void init_k_means(MImage &X_s, float *means, int nbClasses) {
	int k = 0;
	//init means
	for (int i = 0; i < nbClasses; i++) {
		means[i] = 0.0f;
	}

	do {
		int x = (std::rand() % 100 * 1.0f) / 100 * X_s.MXSize();
		int y = (std::rand() % 100 * 1.0f) / 100 * X_s.MYSize();
		float color = X_s.MGetColor(x, y);

		bool equal = false;
		if (k != 0) {
			for (int i = 0; i < nbClasses; i++)
			{
				if (std::abs(means[i] - color) < 1) {
					equal = true;
					break;
				}
			}
		}

		if (!equal) {
			means[k] = X_s.MGetColor(x, y);
			k++;
		}
	} while (nbClasses != k);
}


/*
N-class KMeans segmentation
Resulting values are copied in parameters 'means','stddev', and 'apriori'
The 'apriori' parameter contains the proportion of each class.
The resulting label Field is copied in the current image (this->MImgBuf)
*/
std::vector<std::vector<int> > MImage::MKMeansSegmentation(float *means, float *stddev, float *apriori, int nbClasses)
{
	init_k_means(*this, means, nbClasses);
	std::vector<std::vector<int> > map_black_white(MXS, std::vector<int>(MYS, 0));

	//Calculate means while mean not converge 
	bool same_means = true;
	do {
		same_means = true;
		std::vector<std::vector<float> > mu_classes(nbClasses);
		for (int y = 0; y < MYSize(); y++)
		{
			for (int x = 0; x < MXSize(); x++)
			{
				//find mean of class more near than point (x,y)
				int k = 0;
				for (int i = 1; i < nbClasses; i++) {
					if (abs(means[k] - MImgBuf[x][y].r) > std::abs(means[i] - MImgBuf[x][y].r)) {
						k = i;
					}
				}
				mu_classes[k].push_back(MImgBuf[x][y].r);
				map_black_white[x][y] = k;
			}
		}

		//calculate new means, new stddev and new apriori
		for (int k = 0; k < nbClasses; k++) {
			float new_means = 0.0f;
			stddev[k] = 0.0f;

			//means
			for (int i = 0; i < mu_classes[k].size(); i++) {
				new_means += mu_classes[k][i];
			}
			new_means = new_means / mu_classes[k].size();
			//stddev
			for (int i = 0; i < mu_classes[k].size(); i++) {
				stddev[k] += pow(mu_classes[k][i] - new_means, 2);
			}
			stddev[k] = sqrt(stddev[k] / mu_classes[k].size());
			//apriori
			apriori[k] = (mu_classes[k].size() * 1.0f) / (MXS * MYS);

			if (fabs(means[k] - new_means) > std::numeric_limits<float>::epsilon()) {
				same_means = false;
			}
			means[k] = new_means;
		}
	} while (!same_means);

	//new image
	for (int y = 0; y < MYSize(); y++)
	{
		for (int x = 0; x < MXSize(); x++)
		{
			MImgBuf[x][y].r = ((map_black_white[x][y] * 1.0f) / (nbClasses - 1)) * 255;
		}
	}
	return map_black_white;
}

/*
N-class Soft KMeans segmentation
Resulting values are copied in parameters 'means' and 'stddev'.
The 'apriori' parameter contains the proportion of each class.
The resulting label Field is copied in the current image (this->MImgBuf)
*/
void MImage::MSoftKMeansSegmentation(float *means, float *stddev, float *apriori, float beta, int nbClasses)
{
	init_k_means(*this, means, nbClasses);
	std::vector<std::vector<std::vector<float> > > map_black_white(MXS, std::vector<std::vector<float> >(MYS, std::vector<float>(nbClasses, 0.0f)));

	//Calculate means while mean not converge
	bool same_means = true;
	do {
		same_means = true;
		//Calculate d_c
		for (int y = 0; y < MYSize(); y++)
		{
			for (int x = 0; x < MXSize(); x++)
			{
				for (int k = 0; k < nbClasses; k++) {
					map_black_white[x][y][k] = exp(-(beta * std::abs(MImgBuf[x][y].r - means[k])));
				}
			}
		}

		//finished calculate part of P(c|x_s), sum(e^-d_r)
		std::vector<float> new_means(nbClasses, 0.0f);
		std::vector<float> sum_denum_means(nbClasses, 0.0f);
		for (int y = 0; y < MYSize(); y++)
		{
			for (int x = 0; x < MXSize(); x++)
			{
				//calcul d_r
				float d_r = 0.0f;
				for (int i = 0; i < nbClasses; i++)
					d_r += map_black_white[x][y][i];

				for (int k = 0; k < nbClasses; k++) {
					//P(c|x_s)
					map_black_white[x][y][k] = map_black_white[x][y][k] / d_r;
					//means
					new_means[k] += (map_black_white[x][y][k] * MImgBuf[x][y].r);
					sum_denum_means[k] += map_black_white[x][y][k];
				}
			}
		}
		//u_c
		for (int k = 0; k < nbClasses; k++) {
			new_means[k] = new_means[k] / sum_denum_means[k];

			//init apriori and stddev
			apriori[k] = 0.0f;
			stddev[k] = 0.0f;
		}

		//apriori
		std::vector<std::vector<float> > mu_classes(nbClasses);
		for (int y = 0; y < MYSize(); y++)
		{
			for (int x = 0; x < MXSize(); x++)
			{
				int c = 0;
				for (int k = 1; k < nbClasses; k++) {
					if (fabs(map_black_white[x][y][c] - map_black_white[x][y][k]) < 0.5)
						c = k;
				}
				mu_classes[c].push_back(MImgBuf[x][y].r);
				apriori[c] += 1;
			}
		}

		//stddev and finish apriori
		for (int k = 0; k < nbClasses; k++) {
			for (int i = 0; i < mu_classes[k].size(); i++) {
				//(x-mu)^2
				stddev[k] += pow(mu_classes[k][i] - new_means[k], 2);
			}
			if (mu_classes[k].size() != 0)
				stddev[k] = sqrt(stddev[k] / mu_classes[k].size());
			apriori[k] = apriori[k] / (MXS * MYS);
			if (fabs(means[k] - new_means[k]) > 0.5) {
				//std::cout << means[k] << " ; " << new_means[k] << std::endl;
				same_means = false;
			}
			//std::cout << std::endl;
			means[k] = new_means[k];
		}
	} while (!same_means);
}


/*
N-class Expectation maximization segmentation

init values are in 'means', 'stddev' and 'apriori'
Resulting values are copied in parameters 'means', 'stddev' and 'apriori'
The 'apriori' parameter contains the proportion of each class.
The resulting label Field is copied in the current image (this->MImgBuf)
*/
void MImage::MExpectationMaximization(float *means, float *stddev, float *apriori, int nbClasses)
{
}

void calculate_w(int *w, std::vector<std::vector<int>> &map_black_white, int xSeed, int ySeed, float *means, int nbClasses) {
	std::vector<std::pair<int, int>> x_y;
	x_y.push_back(std::pair<int, int>(xSeed + 1, ySeed));
	x_y.push_back(std::pair<int, int>(xSeed - 1, ySeed));
	x_y.push_back(std::pair<int, int>(xSeed, ySeed + 1));
	x_y.push_back(std::pair<int, int>(xSeed, ySeed - 1));

	x_y.push_back(std::pair<int, int>(xSeed + 1, ySeed + 1));
	x_y.push_back(std::pair<int, int>(xSeed - 1, ySeed - 1));
	x_y.push_back(std::pair<int, int>(xSeed + 1, ySeed - 1));
	x_y.push_back(std::pair<int, int>(xSeed - 1, ySeed + 1));

	for (int k = 0; k < nbClasses; k++) {
		w[k] = 0;
	}

	for (auto i = x_y.begin(); i != x_y.end(); i++) {
		if (i->first >= 0 && i->first < map_black_white.size() && i->second >= 0 && i->second < map_black_white[0].size()) {
			for (int k = 0; k < nbClasses; k++) {
				if (map_black_white[i->first][i->second] != k)
					w[k] += 1;
			}
			//find mean of class more near than point (x,y)
			/*int k = 0;
			for (int j = 1; j < nbClasses; j++) {
				if (abs(means[k] - X_s.MGetColor(i->first, i->second)) > std::abs(means[j] - X_s.MGetColor(i->first, i->second))) {
					k = j;
				}
			}
			w[k] += 1;*/
		}
	}
}

/*
N-class ICM segmentation
beta : Constant multiplying the apriori function
The resulting label Field is copied in the current image (this->MImgBuf)
*/
void MImage::MICMSegmentation(float beta, int nbClasses)
{
	MImage img_prev, img_orig;
	float *means;
	float *stddev;
	float *apriori;
	means = new float[nbClasses];
	stddev = new float[nbClasses];
	apriori = new float[nbClasses];
	img_orig = *this;
	std::vector<std::vector<int>> map_black_white = MKMeansSegmentation(means, stddev, apriori, nbClasses);

	/*img_prev = *this;
	printf(" Means : "); for (int i = 0; i<nbClasses; i++) printf(" %f", means[i]);
	printf("\n Stddev : "); for (int i = 0; i<nbClasses; i++) printf(" %f", stddev[i]);
	printf("\n Prior : "); for (int i = 0; i<nbClasses; i++) printf(" %f", apriori[i]);
	printf("\n");
	img_prev.MRescale();
	img_prev.MSaveImage("outKMeans.pgm", PGM_ASCII);*/
	
	bool same_means = true;
	do{
		//img_prev = *this;
		same_means = true;
		int nb_pxl_diff = 0;
		std::vector<std::vector<float> > mu_classes(nbClasses);
		for (int y = 0; y < MYSize(); y++)
		{
			for (int x = 0; x < MXSize(); x++)
			{
				float *tag;
				tag = new float[nbClasses];

				int tag_min = 0;
				int *w;
				w = new int[nbClasses];
				calculate_w(w, map_black_white, x, y, means, nbClasses);
				for (int k = 0; k < nbClasses; k++) {
					float u = -std::log(exp(-pow(img_orig.MImgBuf[x][y].r * 1.0f - means[k], 2) / (2 * pow(stddev[k], 2))) / (stddev[k] * pow(2 * M_PI, 0.5)));
					tag[k] = u + (w[k] * beta);
					if (tag_min != k && (tag[k] - tag[tag_min]) < 0.0f) {
						tag_min = k;
					}
				}
				//MImgBuf[x][y].r = ((tag_min * 1.0f) / (nbClasses - 1)) * 255.0f;
				mu_classes[tag_min].push_back(img_orig.MImgBuf[x][y].r);
				if (tag_min != map_black_white[x][y]) {
					same_means = false;
					nb_pxl_diff++;
				}
				map_black_white[x][y] = tag_min;
			}
		}
		printf("\n diff pxl: %i", nb_pxl_diff);

		//calculate new means, new stddev and new apriori
		for (int k = 0; k < nbClasses; k++) {
			float new_means = 0.0f;
			stddev[k] = 0.0f;

			//means
			for (int i = 0; i < mu_classes[k].size(); i++) {
				new_means += mu_classes[k][i];
			}
			new_means = new_means / mu_classes[k].size();
			//stddev
			for (int i = 0; i < mu_classes[k].size(); i++) {
				stddev[k] += pow(mu_classes[k][i] - new_means, 2);
			}
			stddev[k] = sqrt(stddev[k] / mu_classes[k].size());
			//apriori
			apriori[k] = (mu_classes[k].size() * 1.0f) / (MXS * MYS);

			if (fabs(means[k] - new_means) > std::numeric_limits<float>::epsilon()) {
				same_means = false;
			}
			means[k] = new_means;
		}

	} while (!same_means);

	//new image
	for (int y = 0; y < MYSize(); y++)
	{
		for (int x = 0; x < MXSize(); x++)
		{
			MImgBuf[x][y].r = ((map_black_white[x][y] * 1.0f) / (nbClasses - 1)) * 255;
		}
	}
}

/*
N-class Simulated annealing segmentation
beta : Constant multiplying the apriori function
Tmax : Initial temperature (initial temperature)
Tmin : Minimal temperature allowed (final temperature)
coolingRate : rate by which the temperature decreases
The label Field copied in the current image (this->MImgBuf)
*/
void MImage::MSASegmentation(float beta, float Tmin, float Tmax, float coolingRate, int nbClasses)
{
}



/*
Interactive graph cut segmentation

the implementation is inspired of the following paper
Y Boykov and M-P Jolly "Interactive Graph Cuts for Optimal Boundary & Region Segmentation of Objects in N-D images".
In International Conference on Computer Vision, (ICCV), vol. I, pp. 105-112, 2001.
The resulting label Field is copied in the current image (this->MImgBuf)
*/
void MImage::MInteractiveGraphCutSegmentation(MImage &mask, float sigma)
{
	//int *result = new int[num_pixels];   // stores result of optimization

	try {
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(MXS, MYS, MZS);

		//// first set up data costs individually
		//for (int i = 0; i < num_pixels; i++)
		//	for (int l = 0; l < num_labels; l++)
		//		if (i < 25) {
		//			if (l == 0) gc->setDataCost(i, l, 0);
		//			else gc->setDataCost(i, l, 10);
		//		}
		//		else {
		//			if (l == 5) gc->setDataCost(i, l, 0);
		//			else gc->setDataCost(i, l, 10);
		//}

		//// next set up smoothness costs individually
		//for (int l1 = 0; l1 < num_labels; l1++)
		//	for (int l2 = 0; l2 < num_labels; l2++) {
		//		int cost = (l1 - l2)*(l1 - l2) <= 4 ? (l1 - l2)*(l1 - l2) : 4;
		//		gc->setSmoothCost(l1, l2, cost);
		//}

		//printf("\nBefore optimization energy is %d", gc->compute_energy());
		//gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		//printf("\nAfter optimization energy is %d", gc->compute_energy());

		//for (int i = 0; i < num_pixels; i++)
		//	result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e) {
		e.Report();
	}

	//delete[] result;
}


/* =================================================================================
====================================================================================
======================              Operators                 ======================
====================================================================================
====================================================================================*/
void MImage::operator= (const MImage &copy)
{
	if (copy.MIsEmpty()) {
		MFreeMemory();
		return;
	}

	if (!MSameSize(copy))
		MAllocMemory(copy.MXS, copy.MYS, copy.MZS);

	for (int y = 0; y<MYS; y++)
		for (int x = 0; x<MXS; x++)
			MImgBuf[x][y] = copy.MImgBuf[x][y];

}

void MImage::operator= (float val)
{
	if (MIsEmpty()) {
		return;
	}

	for (int y = 0; y<MYS; y++)
		for (int x = 0; x<MXS; x++)
			for (int z = 0; z<MZS; z++)
				MSetColor(val, x, y, z);

}

/* =================================================================================
====================================================================================
======================              Ohters                 ======================
====================================================================================
====================================================================================*/


float MImage::MPottsEnergy(const MImage &img, int x, int y, int label) const
{
	return 0.0;
}

float MImage::MComputeGlobalEnergy(const MImage &X, float *mean, float *stddev, float beta, int nbClasses) const
{
	return 0.0;
}