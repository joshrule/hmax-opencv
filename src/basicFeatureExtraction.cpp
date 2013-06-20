/*
 *  basicFeatureExtraction.cpp
 *  hmax_xcode
 *
 *  Created by Drenkow, Nathan G. on 10/23/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

using namespace std;

#ifndef USE_LEGACY_HMAX
#define USE_LEGACY_HMAX 0
#endif

#include <vector>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <dirent.h>

#if 0
#include <sys/time.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#endif

#include "Filters.h"
#include "cHMAX.h"
#include "utilities.h"

#include "opencv2/opencv.hpp"

const double MIN_SIZE = 240;

int main(int argc, char* argv[]) {

	string imageDir = "";
	string patchDir = "";
	string patchFile = "";
	
    	int i=1; //Recall: argv[0] is the executable call, so argv[1] is the first desired command line argument
	while (i<argc) {
		string arg(argv[i]);
		
		if (arg == "-id") { //Image directory
			imageDir.assign(argv[i+1]);
		} else if (arg == "-pd") {
            		patchDir.assign(argv[i+1]);
        	} else if (arg == "-pf") {
            		patchFile.assign(argv[i+1]);
		} else {
			cout << "Invalid argument" << endl;
			exit(-1);
		}        
		i += 2;
    	}
    
	CHMAX hmx(patchDir + patchFile);
	
	vector<string> imgFilenames;
	vector<cv::Mat> imgVector;
	
	size_t period = patchFile.find_last_of(".");
	string patchBase = patchFile.substr(0,int(period));

	getFileListFromDir(imageDir, imgFilenames);
    
	for (int i=0; i<imgFilenames.size(); i++) {
        cv::Mat img;
        cv::Mat src;
        readFileImages(imgFilenames[i], img);
//        readimages_resize(src,MIN_SIZE,img);
        
		cv::Mat features;
		hmx.runHMAX(img, features);
		string featureFile = imgFilenames[i].substr(0,int(imgFilenames[i].find("."))) + ".xml";
		writeFeaturesToFile(featureFile,patchBase,features);
		
//		cv::Mat temp2;
//		readFeaturesFromFile(featureFile, patchBase, temp2);
	}
	
}//end main
