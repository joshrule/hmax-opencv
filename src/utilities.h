#ifndef __UTILITIES__
#define __UTILITIES__

#include "opencv2/opencv.hpp"
// #include <opencv2/gpu/gpu.hpp>

void matlab_rgb2gray(const string& filename,cv::Mat& img);
void matlab_rgb2gray(const cv::Mat& color_img,cv::Mat& img);

//========== IMAGE RESIZING ===========
void matlab_imresize(const cv::Mat& src,const double scale,const int resizeAlongCols,cv::Mat& dst);
void matlab_imresize_full(const cv::Mat& src,const double scale,cv::Mat& dst);

void readimages_resize(const cv::Mat& src,const double minsize,cv::Mat& dst);
void readimages_resize(const vector<cv::Mat> &src_vect,const double minsize, vector<cv::Mat> &dst_vect);

void cubic(const cv::Mat& x, cv::Mat& y, double scale = 1);

//=====================================
void readFileNamesFromText(const string fileName, const string positiveClassName, vector<string> &positiveFileNames, vector<string> &negativeFileNames);
void readFileNamesFromText(const string fileName, vector<string> &imgFilenames);

// Method assumes FileNames vector strings do NOT contain full image paths
void readFileImages(const string imgdir, vector<string> &FileNames, vector<cv::Mat> &Images);

// Method assumes FileNames strings DO contain full image paths and retrieves all specified images
void readFileImages(const vector<string> &FileNames, vector<cv::Mat> &Images);
void readFileImages(const string &FileNames, cv::Mat &Images);//PAR overloaded function to just read one image

// Method assumes fileName contains one full image filepath per line 
void getFileImages(const string fileName, vector<cv::Mat> &Images);

// Method retrieves a list of valid images in specified directory
void getFileListFromDir(const string imgdir, vector<string> &FileNames);

// Method returns a cv::Mat mask that is the same type as the input Mat 
cv::Mat getMask(const cv::Mat& in, double thresh = 0, int compareMethod = cv::CMP_GE);

// Method reads a 2-column file (could be modified later) and returns vectors corresponding to columns
void readTwoColFile(const string filename, vector<string> &col1, vector<string> &col2);

// Method reads a 2-column file (could be modified later) and a vector and Mat corresponding to columns
void readTwoColFile(const string filename, vector<string> &col1, cv::Mat &col2);

// Method writes features to xml file for later use
void writeFeaturesToFile(const string filename, const string source, const cv::Mat& features);

// Method reads features from xml file
void readFeaturesFromFile(const string filename, const string patchType, cv::Mat& features);

// Method reads features from multiple xml files
void readFeaturesFromFile(const vector<string> filename, const string patchType, vector<cv::Mat>& features);

// Method checks whether feature vector for given patch type has been extracted already
bool existFeaturesInXml(const string imgFilename, const string patchType);

// Method adds directory path to filenames
vector<string> fullfile(const string dir, const vector<string> filenames);
string fullfile(const string dir, string filename);

#endif
