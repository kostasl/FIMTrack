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

#ifndef INVERTDIALOG_HPP
#define INVERTDIALOG_HPP

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-register"
#include <QDialog>
#pragma clang diagnostic pop

namespace Ui {
    class InvertDialog;
}

class InvertDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit InvertDialog(QWidget *parent = 0);
    ~InvertDialog();
    
public slots:
    void setBeforeTimeSteps(QStringList const& timeSteps);
    void setAfterTimeSteps(QStringList const& timeSteps);
    
private slots:
    void setupConnections();
    
    void on_pbtOK_clicked();
    
    void on_pbtCancel_clicked();
    
signals:
    void sendParameter(bool before, uint bTime, bool current, bool after, uint aTime);
    
private:
    Ui::InvertDialog *ui;
};

#endif // INVERTDIALOG_HPP
