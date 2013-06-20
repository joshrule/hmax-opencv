#ifndef __MYTIMER__
#define __MYTIMER__

// Classes for MyMouse and MyTimer

//===========================================================================================
//===========================================================================================

// Copy of cv_mouse from cv_utilities
class MyMouse
{
public:
    static void start(const cv::string& a_img_name)
    {
        cvSetMouseCallback(a_img_name.c_str(), MyMouse::cv_on_mouse, 0);
    }
    static int event(void)
    {
        int l_event = m_event;
        m_event = -1;
        return l_event;
    }
    static int x(void)
    {
        return m_x;
    }
    static int y(void)
    {
        return m_y;
    }
    
private:
    static void cv_on_mouse(int a_event, int a_x, int a_y, int, void *)
    {
        m_event = a_event;
        m_x = a_x;
        m_y = a_y;
    }
    
    static int m_event;
    static int m_x;
    static int m_y;
};

//===========================================================================================
//===========================================================================================

// Adapted from cv_timer in cv_utilities
class MyTimer
{
public:
    MyTimer() : start_(0), time_(0) {}
    
    void start()
    {
        start_ = cv::getTickCount();
    }
    
    void stop()
    {
        CV_Assert(start_ != 0);
        int64 end = cv::getTickCount();
        time_ += end - start_;
        start_ = 0;
    }
    
    double time()
    {
        double ret = time_ / cv::getTickFrequency();
        time_ = 0;
        return ret;
    }
    
private:
    int64 start_, time_;
};

//===========================================================================================

#endif