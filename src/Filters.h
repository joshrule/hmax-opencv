#ifndef __FILTERS__
#define __FILTERS__

#include "opencv2/opencv.hpp"
// #include <opencv2/gpu/gpu.hpp>


void conv2(cv::Mat* src, cv::Mat* kernel, cv::Mat* dst);

void sumfilter(cv::Mat* I, cv::Scalar radius, cv::Mat* dst);

void WindowedPatchDistance(cv::Mat& Im, cv::Mat& Patch, cv::Mat& dst);

void WindowedPatchDistance(vector<cv::Mat>& Im, vector<cv::Mat>& Patch, cv::Mat& dst);//using vector version

void WindowedPatchDistance(cv::Mat& Im, cv::Mat& Patch, cv::Mat& dst, int doSkipSquareRoot); // version with doSkipSquareRoot flag
void WindowedPatchDistance(vector<cv::Mat>& Im, vector<cv::Mat>& Patch, cv::Mat& dst, int doSkipSquareRoot);//using vector version
void mymaxfilter(cv::Mat& Im, cv::Scalar poolrange, cv::Mat& dst);


#endif
