/*  utilities.cpp: contains varius useful functionality
 *  readFileNamesAndClassesFromText.cpp
 *  FeaturePtTester
 *
 *  Created by Pedro A. Rodriguez
 *  Copyright 2012 Johns Hopkins University APL. All rights reserved.
 *
 *
 */

#include <vector>

using namespace std;

#include <iostream>
#include <fstream>
#include <dirent.h>

#if 0
#include <string>
#include <sstream>
#endif

#include "opencv2/opencv.hpp"

#include "utilities.h"

#define DEBUG_FLAG

//=====================================================================================================
// Make two vectors of file names from a fold text file, one of positive file names, one of negative file names
// In the text file, in each line there is an image name followed by the class name for positive images
// or followed by the word negative or another word for negative images
void matlab_rgb2gray(const string& filename,cv::Mat& img) {
	
	cv::Mat img2 = cv::imread(filename); //Read as brg
										 //cv::cvtColor(img2,img2,CV_BGR2GRAY);  //OpenCV function uses different scale factors to go to gray 0.299,0.5870,0.1140 than Matlab
	if (img2.empty()) {
		cout << "Image with file name" << filename << " cannot be opened" << endl; // dfd: return empty-handed instead of crashing
		return;
	}
	img2.convertTo(img2,CV_64FC3);  //Convert to double to do calculations
	cv::Mat img_gray(img2.rows,img2.cols,CV_64FC1,0.0);  //Define new gray image
	
	//Split rgb channels into std vector
	vector<cv::Mat> img_vect;
	cv::split(img2,img_vect);
	
	//Do Weighted sum to go from bgr to gray:
	//img_gray=img_vect[0]*0.1140 + img_vect[1]*0.5870 + img_vect[2]*0.2989; //Matlab help
	//img_gray=img_vect[0]*0.1140 + img_vect[1]*0.5870 + img_vect[2]*0.299; //OpenCV
	img_gray=img_vect[0]*0.114020904255103 + img_vect[1]*0.587043074451121 + img_vect[2]*0.298936021293776; //Actual Matlab
    //OR:
    //cv::addWeighted(img_vect[0],0.1140,img_vect[1],0.5870,0,img_gray);
    //cv::addWeighted(img_gray,1.0,img_vect[2],0.2989,0,img_gray);
    
	// Round using for loop or CV_8UC1 conversion:
	
	//img_gray.convertTo(img_gray,CV_8UC1);
	//OR:
	
    for (int i_chip = 0; i_chip < img_gray.rows; i_chip++) {
		for (int j_chip = 0; j_chip < img_gray.cols; j_chip++) {    
			img_gray.at<double>(i_chip,j_chip)=(double)round(img_gray.at<double>(i_chip,j_chip));
			//cout << "round: " << img_gray.at<double>(i_chip,j_chip) << " " << endl;
		}
    }
	
	//Save final image and divide by 255
	//	img_gray.copyTo(img); // dfd: not needed
	//img.convertTo(img,CV_64F);
	//	cv::divide(img,255.0,img);
	cv::divide(img_gray,255.0,img);
	
	//	img2.release();
	//	img_gray.release();
	
}

//=====================================================================================================
//This version inputs the actual color image:
void matlab_rgb2gray(const cv::Mat& color_img,cv::Mat& img) {
	
	//cv::Mat img2 = cv::imread(filename); //Read as brg
	//cv::cvtColor(img2,img,CV_BGR2GRAY);  //OpenCV function uses different scale factors to go to gray 0.299,0.5870,0.1140 than Matlab
	
	cv::Mat img2; 
	color_img.copyTo(img2);
	img2.convertTo(img2,CV_64FC3);  //Convert to double to do calculations
	cv::Mat img_gray(img2.rows,img2.cols,CV_64FC1,0.0);  //Define new gray image
	
	//Split rgb channels into std vector
	vector<cv::Mat> img_vect;
	cv::split(img2,img_vect);
	
	//Do Weighted sum to go from bgr to gray:
	//img_gray=img_vect[0]*0.1140 + img_vect[1]*0.5870 + img_vect[2]*0.2989; //Matlab help
	//img_gray=img_vect[0]*0.1140 + img_vect[1]*0.5870 + img_vect[2]*0.299; //OpenCV
	img_gray=img_vect[0]*0.114020904255103 + img_vect[1]*0.587043074451121 + img_vect[2]*0.298936021293776; //Actual Matlab    
    //OR:
    //cv::addWeighted(img_vect[0],0.1140,img_vect[1],0.5870,0,img_gray);
    //cv::addWeighted(img_gray,1.0,img_vect[2],0.2989,0,img_gray);
    
	// Round using for loop or CV_8UC1 conversion:
	
	//img_gray.convertTo(img_gray,CV_8UC1);
	//OR:
	
    for (int i_chip = 0; i_chip < img_gray.rows; i_chip++) {
		for (int j_chip = 0; j_chip < img_gray.cols; j_chip++) {    
			img_gray.at<double>(i_chip,j_chip)=(double)round(img_gray.at<double>(i_chip,j_chip));
			//cout << "round: " << img_gray.at<double>(i_chip,j_chip) << " " << endl;
		}
    }
	
	//Save final image and divide by 255
	//	img_gray.copyTo(img); // dfd: not needed
	//img.convertTo(img,CV_64F);
	//	cv::divide(img,255.0,img);
	cv::divide(img_gray,255.0,img);
	
	//	img2.release(); // dfd: not needed
	//	img_gray.release();
	
}


//======================================================================================================
//												IMAGE RESIZING
//======================================================================================================
void matlab_imresize(const cv::Mat& src,const double scale,const int resizeAlongCols,cv::Mat& dst) {
	// NOTE: Height is synonymous with rows, width is synonymous with cols
	double kernel_width = 4;
	
	cv::Size orig_sizes,new_sizes;
	
	cv::Mat img = src.clone();
	img.convertTo(img, CV_64F);
	
	img = (resizeAlongCols) ? img.t() : img;
	
	orig_sizes = img.size();
	
	double width = (double) orig_sizes.width;
	double height = (double) orig_sizes.height;
	
	// Only scale along the specified dimension
	int new_width = (int) width;
	int new_height = (int)ceil(height*scale);
	
	new_sizes.width = new_width;
	new_sizes.height = new_height;
	
	int out_length = new_sizes.height;
	
	kernel_width = (scale < 1) ? kernel_width/scale : kernel_width; //This matches Matlab's scaling for antialiasing when scale < 1 
	
	cv::Mat x(out_length,1,CV_64F);
	for (int i=0; i<out_length; i++) {
		x.at<double>(i,0) = i+1;
	}
	
	cv::Mat u(x.rows, x.cols, CV_64F);
	u = (x/scale) + 0.5 * (1 - (1/scale));
	
	cv::Mat left(u.rows, u.cols, CV_64F);
	left = u - kernel_width/2;
	for (int i=0; i<left.rows; i++) {
		left.at<double>(i,0) = floor(left.at<double>(i,0));
	}
	
	double P = ceil(kernel_width) + 2;
	
	cv::Mat rowP(1,(int)P,CV_64F);
	for (int i=0; i<(int)P; i++) {
		rowP.at<double>(0,i) = i;
	}
	
	// Compute weights and indices
	cv::Mat indices(out_length,(int)P,CV_64F);
	cv::Mat indicesUnclamped(out_length,(int)P,CV_64F);
	cv::Mat weights;
	for (int i=0; i<out_length; i++) {
		for (int j=0; j<(int)P; j++) {
			indices.at<double>(i,j) = left.at<double>(i,0) + rowP.at<double>(0,j);  
			indicesUnclamped.at<double>(i,j) = indices.at<double>(i,j); //Includes negatives and values greater than the original length
			indices.at<double>(i,j) = max(indices.at<double>(i,j), 1.0);  //Remove negative indices
			indices.at<double>(i,j) = min(indices.at<double>(i,j),height); //Remove indices greater than original length
		}		
		cv::Mat tempIdx(1,indicesUnclamped.cols,CV_64F);
		tempIdx = indicesUnclamped.row(i).clone();
		cv::Mat tempU(1,u.cols,CV_64F);
		tempU = u.row(i).clone();
		cv::Mat tempWts;
		cubic(tempU - tempIdx,tempWts,scale);
		
		cv::Scalar wtSum = sum(tempWts);
		tempWts = tempWts / ((double) wtSum.val[0]);
		weights.push_back(tempWts);
	}
	
	//Account for zero-indexing
	indicesUnclamped = indicesUnclamped - 1;
	indices = indices - 1;
	
	dst.create(new_sizes.height,new_sizes.width,CV_64F);
	
	// Do actual resizing
	for (int iCol=0; iCol<new_sizes.width; iCol++) {
		
		for (int iRow=0; iRow < out_length; iRow++){
			double pixelval = 0;
			for (int j=0; j<indices.cols; j++) {
				int oldRow = (int) indices.at<double>(iRow,j);
				pixelval += weights.at<double>(iRow,j)*img.at<double>(oldRow,iCol);
			}
			dst.at<double>(iRow,iCol) = pixelval;
		}
	}	
	
	dst = (resizeAlongCols) ? dst.t() : dst;	
}

//=====================================================================================================
void matlab_imresize_full(const cv::Mat& src,const double scale,cv::Mat& dst) {
	//Rescales images using matlab_imresize to scale columns and rows.
	
	cv::Size orig_sizes;
	orig_sizes = src.size();
	
	double width=(double)orig_sizes.width;
	double height=(double)orig_sizes.height;
	
	double min_side = min(width,height);
	
	cv::Mat temp;
	if (min_side == width) {
		matlab_imresize(src,scale,1,temp); //Resize along columns first
		matlab_imresize(temp,scale,0,dst);
	} else {
		matlab_imresize(src,scale,0,temp);
		matlab_imresize(temp,scale,1,dst); //Resize along columns second
	}	
	
}

//=====================================================================================================
void readimages_resize(const cv::Mat& src,const double minsize,cv::Mat& dst) {
	
	cv::Size orig_sizes;
	orig_sizes = src.size();
	
	double width=(double)orig_sizes.width;
	double height=(double)orig_sizes.height;
	
	double min_side = min(width,height);
	
	if (min_side > minsize) {
		double scale = minsize/min_side;
		cv::Mat temp;
		if (min_side == width) {
			matlab_imresize(src,scale,1,temp); //Resize along columns first
			matlab_imresize(temp,scale,0,dst);
		} else {
			matlab_imresize(src,scale,0,temp);
			matlab_imresize(temp,scale,1,dst); //Resize along columns second
		}		
		
	} else {
		dst = src;
	}
	
}

//=====================================================================================================
void readimages_resize(const vector<cv::Mat> &src_vect,const double minsize, vector<cv::Mat> &dst_vect) {
	
	double width,length,min_side;
	cv::Mat src,dst;
	cv::Size orig_sizes;
	
	int numImgs = src_vect.size();
	
	for (int i = 0; i < numImgs; i++) {
		
		src_vect[i].copyTo(src);
		
		orig_sizes=src.size();
		
		width=(double)orig_sizes.width;
		length=(double)orig_sizes.height;
		
		min_side=min(width,length);
		
		if (min_side > minsize) {
			double scale = minsize/min_side;
			cv::Mat temp;
			if (min_side == width) {
				matlab_imresize(src,scale,1,temp); //Columns first
				matlab_imresize(temp,scale,0,dst);
			} else {
				matlab_imresize(src,scale,0,temp);
				matlab_imresize(temp,scale,1,dst); //Columns second
			}		
		} else {
			dst = src;
		}
		
		dst_vect.push_back(dst);	
	}
	
}

//=====================================================================================================
// Interpolation kernel without antialiasing
void cubic(const cv::Mat& x_in, cv::Mat& y, double scale){
	cv::Mat x = x_in.clone();
	x.convertTo(x,CV_64F);
	
	x = (scale < 1) ? scale*x : x;  //This matches Matlab's scaling for antialiasing when scale < 1
	
	cv::Mat absx = abs(x);
	cv::Mat absx2 = absx.mul(absx);
	cv::Mat absx3 = absx2.mul(absx);
	
	cv::Mat term1 = (1.5*absx3 - 2.5*absx2 + 1);
	cv::Mat term1idx = getMask(absx,1,cv::CMP_LE);
	term1 = term1idx.mul(term1);
	
	cv::Mat term2 = (-0.5*absx3 + 2.5*absx2 - 4*absx + 2);
	cv::Mat term2idx;
	multiply(getMask(absx,1,cv::CMP_GT),getMask(absx,2,cv::CMP_LE),term2idx);
	term2 = term2idx.mul(term2);
	
	y = (scale < 1) ? scale*(term1 + term2) : (term1 + term2); //This matches Matlab's scaling for antialiasing when scale < 1
	y.convertTo(y,CV_64F);
}

//======================================================================================================
//										FILE I/O
//======================================================================================================
void readFileNamesFromText(const string fileName, const string positiveClassName, vector<string> &positiveFileNames, vector<string> &negativeFileNames)
{
	ifstream ifs;
	std::string str;
	std::string imageName, className;
	
	ifs.open(fileName.c_str());
    if (!ifs){
        cout << "File passed to readFileNamesFromText does not exist." << endl;
        exit(1);
    }
	while(!ifs.eof()){
		getline(ifs, str); // also try ifs >> imageName >> className;
        if (str.size() < 3){ // dfd: there should be at least 3 characters such as .jpg
            continue;  // dfd: discard blank lines
        }
		istringstream lineStream(str); // stream
		lineStream >> imageName;
		lineStream >> className;
		if (className.compare(positiveClassName)==0){ // image is positive
			positiveFileNames.push_back(imageName);
		}
		else {
			negativeFileNames.push_back(imageName);
		}
	}
	ifs.close();
}

//=====================================================================================================
void readFileNamesFromText(const string fileName, vector<string> &imgFilenames)
{
	ifstream ifs;
	std::string str;
	std::string imageName, className;
	ifs.open(fileName.c_str());
    if (!ifs){
        cout << "File passed to readFileNamesFromText does not exist." << endl;
        exit(1);
    }
	while(!ifs.eof()){
		getline(ifs, str); // also try ifs >> imageName >> className;
        if (str.size() < 3){ // dfd: there should be at least 3 characters such as .jpg
            continue;  // dfd: discard blank lines
        }
		istringstream lineStream(str); // stream
		lineStream >> imageName;
		imgFilenames.push_back(imageName);
	}
	ifs.close();
}

//=====================================================================================================
void readFileImages(const string imgdir, vector<string> &FileNames, vector<cv::Mat> &Images)
{
	//Inputs:
	// image directory= where images are located
	// positiveFileNames = vector of strings with the names of the images
	// Output:
	// Images = vector of grayscale images
	
	int numImgs = FileNames.size();
	
	for (int i = 0; i < numImgs; i++) {
		cv::Mat temp_img;
		string fname = imgdir + FileNames[i];
		matlab_rgb2gray(fname,temp_img);
		if (temp_img.empty()){ // dfd: if image could not be found in directory, keep going instead of crashing
			continue;
		}
		Images.push_back(temp_img);
		
	}	
}

//=====================================================================================================
void readFileImages(const vector<string> &FileNames, vector<cv::Mat> &Images)
{
	//Inputs:
	// image directory= where images are located
	// positiveFileNames = vector of strings with the names of the images
	// Output:
	// Images = vector of grayscale images
	
	int numImgs = FileNames.size();
	
	for (int i = 0; i < numImgs; i++) {
		cv::Mat temp_img;
		string fname = FileNames[i];
		matlab_rgb2gray(fname,temp_img);
		if (temp_img.empty()){ // // dfd: if image could not be found in directory, keep going instead of crashing
			continue;
		}
		Images.push_back(temp_img);
		
	}
	
}

//PAR Overloaded function to do read one image at the time:
void readFileImages(const string &FileNames, cv::Mat &Images)
{
	//Inputs:
	// image directory= where images are located
	// positiveFileNames = vector of strings with the names of the images
	// Output:
	// Images = vector of grayscale images
	
	int numImgs = 1;//FileNames.size();
	
	for (int i = 0; i < numImgs; i++) {
		cv::Mat temp_img;
		string fname = FileNames;
		matlab_rgb2gray(fname,temp_img);
		if (temp_img.empty()){ // // dfd: if image could not be found in directory, keep going instead of crashing
			continue;
		}
		Images=temp_img;
		
	}
	
}

//=====================================================================================================
void getFileImages(const string fileName, vector<cv::Mat> &Images) {
	
	vector<string> fileNames;
	readFileNamesFromText(fileName, fileNames);
	
	readFileImages(fileNames,Images);	
}

//=====================================================================================================
//					Extract Image List from Directory
//=====================================================================================================
void getFileListFromDir(const string imgDir, vector<string> &FileNames) {
	
	char *directory = (char *)imgDir.c_str();
	
	vector<string> validExt;
	validExt.push_back(".png"); validExt.push_back(".jpg");
	validExt.push_back(".tiff"); validExt.push_back(".bmp");
	validExt.push_back(".tif");	validExt.push_back(".jpeg");
	validExt.push_back(".jpe");	validExt.push_back(".dib");
	validExt.push_back(".jp2");	validExt.push_back(".pbm");
	validExt.push_back(".pgm");	validExt.push_back(".ppm");
	validExt.push_back(".sr");	validExt.push_back(".ras");
	
	DIR *dir = opendir(directory);
	if(dir)
	{
		struct dirent *ent;
		while((ent = readdir(dir)) != NULL)
		{
			string filename(ent->d_name);
			cout << filename << endl;
			
			if (filename != "." && filename != "..") {
				for (int i=0; i<validExt.size(); i++) {
					// Check if current file is an image
					if (filename.find(validExt[i]) != string::npos) {
						FileNames.push_back(imgDir + filename);
					}
				}
			}
		}		
	} 
	else {
		fprintf(stderr, "Error opening directory\n");
	}
	
	
} 

//==============================================================================================================================
cv::Mat getMask(const cv::Mat& in, double thresh, int compareMethod) {
	
	cv::Mat tempmask(in.rows,in.cols,CV_8UC1);
	
	cv::compare(in,thresh,tempmask,compareMethod);
	cv::divide(tempmask,255.0,tempmask);
	tempmask.convertTo(tempmask,in.type());
	return tempmask;
}


//=====================================================================================================
void writeFeaturesToFile(const string filename, const string source, const cv::Mat& features) {
	cv::FileStorage fs(filename, cv::FileStorage::READ);
	
	cv::FileNode sourceNode = fs[source];
	if (sourceNode.empty()) {
		fs.release();
		fs.open(filename, cv::FileStorage::APPEND);
		fs << source << features;		
	} else {
		//Do nothing.  Features for provided field name already exist.
	}
	
	fs.release();
	
}

//=====================================================================================================
void readTwoColFile(const string filename, vector<string> &col1, vector<string> &col2) {
	ifstream ifs;
	std::string str;
	std::string col1str, col2str;
	
	ifs.open(filename.c_str());
    if (!ifs){
        cout << "File passed to readTwoColFile does not exist." << endl;
        exit(1);
	}
	while(!ifs.eof()){
		getline(ifs, str); 
		istringstream lineStream(str); // stream
		getline(lineStream,col1str,'\t');
		getline(lineStream,col2str,'\t');
		col1.push_back(col1str);
		col2.push_back(col2str);
	}
	ifs.close();
}

//=====================================================================================================
void readTwoColFile(const string filename, vector<string> &col1, cv::Mat &col2) {
	ifstream ifs;
	std::string str;
	std::string col1str, col2str;
	
	ifs.open(filename.c_str());
    if (!ifs){
        cout << "File passed to readTwoColFile does not exist." << endl;
        exit(1);
	}
	
	while(!ifs.eof()){
		getline(ifs, str); 
		istringstream lineStream(str); // stream
		lineStream >> col1str;
		lineStream >> col2str;
		
		col1.push_back(col1str);
		
		double col2val;
		stringstream ss(col2str); //turn the string into a stream
		ss >> col2val; //convert
		
		col2.push_back(col2val);
	}
	ifs.close();
}

//=====================================================================================================
// NOTE: Method replaces filename extension with .xml
void readFeaturesFromFile(const string filename, const string patchType, cv::Mat& features) {
	size_t period = filename.find_last_of(".");
	string filebase = filename.substr(0, int(period));
	cv::FileStorage fs((filebase + ".xml"), cv::FileStorage::READ);

	if (fs.isOpened()) {
		fs[patchType] >> features;
	}	
}

//=====================================================================================================
// NOTE: Method replaces filename extension with .xml
void readFeaturesFromFile(const vector<string> filenames, const string patchType, vector<cv::Mat>& features){
	cv::FileStorage fs;
	
	for (int i=0; i<filenames.size(); i++) {
		size_t period = filenames[i].find_last_of(".");
		string filebase = filenames[i].substr(0, int(period));

		fs.open((filebase + ".xml"), cv::FileStorage::READ);
		cv::Mat temp;
		if (fs.isOpened()) {
			fs[patchType] >> temp;		
		}
		features.push_back(temp);
		
		fs.release();
	}
}

//=====================================================================================================
bool existFeaturesInXml(const string imgFilename, const string patchType){
	size_t period = imgFilename.find_last_of(".");
	string filebase = imgFilename.substr(0, int(period));
	cv::FileStorage fs((filebase + ".xml"), cv::FileStorage::READ);
	
	cv::FileNode patchNode = fs[patchType];
	return !patchNode.empty();
}

//=====================================================================================================
vector<string> fullfile(const string dir, const vector<string> filenames){
	vector<string> fullfiles;
	for (int i=0; i<filenames.size(); i++) {
		fullfiles.push_back(dir + filenames[i]);
	}
	return fullfiles;
}

//=====================================================================================================
string fullfile(const string dir, const string filename){
	return (dir + filename);
}