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

#ifndef ARROW_HPP
#define ARROW_HPP

#include "Configuration/FIMTrack.hpp"
#include <cmath>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-register"
#include <QGraphicsLineItem>
#include <QPainter>
#pragma clang diagnostic pop


class Arrow : public QGraphicsLineItem
{
private:
    QColor      mColor;
    
    /**
     * @brief getAngle Obtain the angle in radians between two deltas
     * 12 o'clock is 0.0 * pi 
     * 3 o'clock is 0.5 * pi 
     * 6 o'clock is 1.0 * pi 
     * 9 o'clock is 1.5 * pi
     * @param dx
     * @param dy
     * @return 
     */
    double getAngle(const QLineF &l);
    
public:
    Arrow(QGraphicsItem *parent = 0, QGraphicsScene *scene = 0);
    ~Arrow();
    
    void setColor(QColor const& color) {this->mColor = color;}
    
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    
};

#endif // ARROW_HPP
