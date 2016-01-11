/*****************************************************************************
 * Copyright (c) 2011-2014 The FIMTrack Team as listed in CREDITS.txt        *
 * http://fim.uni-muenster.de                                             	 *
 *                                                                           *
 * This file is part of FIMTrack.                                            *
 * FIMTrack is available under multiple licenses.                            *
 * The different licenses are subject to terms and condition as provided     *
 * in the files specifying the license. See "LICENSE.txt" for details        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FIMTrack is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation, either version 3 of the License, or         *
 * (at your option) any later version. See "LICENSE-gpl.txt" for details.    *
 *                                                                           *
 * FIMTrack is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 * GNU General Public License for more details.                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * For non-commercial academic use see the license specified in the file     *
 * "LICENSE-academic.txt".                                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * If you are interested in other licensing models, including a commercial-  *
 * license, please contact the author at fim@uni-muenster.de      			 *
 *                                                                           *
 *****************************************************************************/

#include "Backgroundsubtractor.hpp"

Backgroundsubtractor::Backgroundsubtractor(t_filePaths const& imagePaths, Undistorter const & undist, QObject *parent) : QObject(parent)
{
    try
    {
        connect(this, SIGNAL(sendLogMessage(QString,LOGLEVEL)), Logger::getInstance(), SLOT(handleLogMessage(QString,LOGLEVEL)));

        if(BackgroundSubstraction::bUseDefault)
        {
            this->iStepSize = 100;
            this->generateBackgroundImage(imagePaths,undist);
        }
        else
        {
            this->generateBackgroundImage(imagePaths, BackgroundSubstraction::iFromImage, BackgroundSubstraction::iOffset, BackgroundSubstraction::iToImage,undist);
        }
    }
    catch (...)
    {
        Logger::getInstance()->addLogMessage(QString("Fatal Error in Backgroundsubtractor::Backgroundsubtractor"), FATAL);
//        Logger::getInstance()->saveLog();
    }
}

cv::Mat& Backgroundsubtractor::substract(cv::Mat const& src, cv::Mat& dst) const
{
    try 
    {
        dst = src - this->backgroundImage;
        return dst;
    }
    catch (...) 
    {
//        emit sendLogMessage(QString("Error in Backgroundsubtractor::substract\nCannot substract backgroundimage"), ERROR);
        return dst;
    }
}

cv::Mat& Backgroundsubtractor::subtractViaThresh(cv::Mat const & src, int const gThresh, cv::Mat & dst) const
{
    try
    {
        /*
         * Alternative implementation of subtract:
         * Since the background image is B(x,y) = min_{t=1,...,T} I(x,y,t), all pixel in
         * src get a bit darker (src(x,y,t) - B(x,y) < src(x,y,t)). Thus the threshold is
         * applied to the result of the subtraction and all pixel above the threshold are
         * used for further calculations:
         * dst(x,y,t) = src(x,y,t)  if src(x,y,t) - B(x,y) > gThresh     or
         *              0           otherwise
         */
//        cv::Mat mask;
//        cv::subtract(src,this->backgroundImage,mask);
//        cv::threshold(mask,mask,gThresh,1,cv::THRESH_BINARY);

//        src.copyTo(dst,mask);
//        return dst;
        src.copyTo(dst);
        for(int row = 0; row < src.rows; ++row) 
        {
            for(int col = 0; col < src.cols; ++col) 
            {
                if(abs(src.at<uchar>(row, col) - this->backgroundImage.at<uchar>(row, col) > abs(gThresh - this->backgroundImage.at<uchar>(row, col))))
                {
                    dst.at<uchar>(row, col) = src.at<uchar>(row, col);
                }
                else
                {
                    dst.at<uchar>(row, col) = 0;
                }
            }
        }
        return dst;
    }
    catch (...)
    {
//        emit sendLogMessage(QString("Error in Backgroundsubtractor::substract\nCannot substract backgroundimage"), ERROR);
        return dst;
    }
}

cv::Mat& Backgroundsubtractor::readImage(std::string const& path, 
                                         cv::Mat& dst)
{
    try 
    {
        dst = cv::imread(path, CV_LOAD_IMAGE_GRAYSCALE);
        return dst;
    } 
    catch (...) 
    {
        emit sendLogMessage(QString("Error in Backgroundsubtractor::readImage\nCould not read image with path: ").append(QString::fromStdString(path)), ERROR);
        dst = cv::Mat();
        return dst;
    }
}

void Backgroundsubtractor::updateBackgroundImage(cv::Mat const& grayImage)
{
    try 
    {
        cv::min(this->backgroundImage, grayImage, this->backgroundImage);
    } 
    catch (...) 
    {
        emit sendLogMessage(QString("Error in Backgroundsubtractor::updateBackgroundImage\nCould update backgroundimage"), ERROR);
    }
}

void Backgroundsubtractor::generateBackgroundImage(std::vector<cv::String> const& imagePaths, Undistorter const & undist)
{
    try
    {
        /* are there some images in the imagelist */
        assert(imagePaths.size() > 0);

        cv::Mat grayImage;
        int iOffset;

        /* calculeate the offset depending on the stepsize */
        double dRatio = (imagePaths.size() / this->iStepSize);
        if(dRatio < 1.0){
            iOffset = 1;
        }else{
            iOffset = static_cast<int>(dRatio);
        }

        /* read the first image as the first backgroundimage */
        this->readImage(imagePaths.at(0), this->backgroundImage);

        if(!this->backgroundImage.empty())
        {
            emit sendLogMessage(QString("Start calculating Backgroundimage with Default Parameters"), INFO);
            emit sendLogMessage(QString("FromImage = ").append(QString::number(0)).append(" Offset = ").append(QString::number(iOffset)).append(" ToImage = ").append(QString::number(imagePaths.size())), DEBUG);
            for (size_t i = 0; i < imagePaths.size(); i += iOffset)
            {
                /* read the next image */
                this->readImage(imagePaths.at(i), grayImage);

                /* if the image could be loaded */
                if(!grayImage.empty())
                {
                    /* update the backgroundimage */
                    this->updateBackgroundImage(grayImage);
                }
            }

            if(undist.isReady())
            {
                cv::Mat dst;
                undist.getUndistortImage(this->backgroundImage,dst);
                dst.copyTo(this->backgroundImage);
            }

            emit sendLogMessage(QString("Calculation of the Backgroundimage Done"), INFO);
        }
    }
    catch (...)
    {
        Logger::getInstance()->addLogMessage(QString("Fatal Error in Backgroundsubtractor::generateBackgroundImage"), FATAL);
//        Logger::getInstance()->saveLog();
    }
}

void Backgroundsubtractor::generateBackgroundImage(const std::vector<cv::String> &imagePaths, unsigned int iFrom, unsigned int iOffset, unsigned int iTo, Undistorter const & undist)
{
    try
    {
        assert(iFrom < imagePaths.size());
        assert(imagePaths.size() > 0);

        cv::Mat grayImage;

        this->readImage(imagePaths.at(0), this->backgroundImage);

        if(!this->backgroundImage.empty())
        {
            emit sendLogMessage(QString("Start calculating Backgroundimage with user Parameters"), INFO);
            emit sendLogMessage(QString("FromImage = ").append(QString::number(iFrom)).append(" Offset = ").append(QString::number(iOffset)).append(" ToImage = ").append(QString::number(iTo)), DEBUG);
            for (size_t i = iFrom; i < iTo; i += iOffset)
            {
                this->readImage(imagePaths.at(i), grayImage);
                if(!grayImage.empty())
                {
                    this->updateBackgroundImage(grayImage);
                }
            }

            if(undist.isReady())
            {
                cv::Mat dst;
                undist.getUndistortImage(this->backgroundImage,dst);
                dst.copyTo(this->backgroundImage);
            }

            emit sendLogMessage(QString("Calculation of the Backgroundimage Done"), INFO);
        }
    }
    catch (...)
    {
        Logger::getInstance()->addLogMessage(QString("Fatal Error in Backgroundsubtractor::generateBackgroundImage"), FATAL);
//        Logger::getInstance()->saveLog();
    }
}
