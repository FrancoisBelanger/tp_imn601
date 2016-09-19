#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "MImage.h"

/*
	"C:\Users\Genevieve\Google Drive\UdS\4e année\6e session\IMN601\TP\prog\satelliteIles.pgm" 90 120 35

	Magic wang image segmentation	
*/ 
int main1(int argc, char **argv)
{
	MImage img;
	
	if(argc<3){
		printf("\n Usage : tp1A image seedX seedY tolerance\n\n"); 
		return 1;
	}

	img.MLoadImage(argv[1]);
	img.MMagicWand(atoi(argv[2]),atoi(argv[3]),atof(argv[4]));
	img.MSaveImage("outMagicWang.pgm",PGM_ASCII);
	
	return 0;
}
