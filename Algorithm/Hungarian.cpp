#include "Hungarian.hpp"

namespace Algorithms
{
    Hungarian::Hungarian()
    {
        this->mIsInit = false;
    }
    
    Hungarian::Hungarian(cv::Mat_<double> const& inputMatrix, 
                         MODE mode)
    {
        this->solve(inputMatrix, mode);
    }
    
    void Hungarian::solve(cv::Mat_<double> const& inputMatrix, 
                          MODE mode)
    {
        this->init(inputMatrix, mode);
        this->runMunkres();
    }
    
    double Hungarian::cost() const
    {
        return this->mCost;
    }  
    
    void Hungarian::init(cv::Mat_<double> const& inputMatrix, 
                         MODE mode)
    {
        this->mIsInit = true;
        this->mCost     = 0.0;
        this->mMaxCost  = 0.0;
        this->mCols     = std::max(inputMatrix.cols, inputMatrix.rows);//inputMatrix.cols;
        this->mRows     = this->mCols;//inputMatrix.rows;
        
        this->mCostMatrixOriginal = inputMatrix;
        
        this->mRowsCovered.clear();
        this->mRowsCovered.resize(this->mRows, false);
        
        this->mColsCovered.clear();
        this->mColsCovered.resize(this->mCols, false);
        
        this->mStarsInCol.clear();
        this->mStarsInCol.resize(this->mCols, 0);
        
        this->mPath.clear();
        this->mPath.resize(std::max(this->mCols, this->mRows) * 2, std::pair<int,int>(-1,-1));
        
        /*
         * Generate the expanded cost matrix by adding extra 0-valued elements in
         * order to make a square matrix. At the same time, find the greatest cost
         * in the matrix (used later if we want to maximize rather than minimize the
         * overall cost.)
         */
//        this->mCostMatrix = cv::Mat::zeros(this->mRows, this->mCols, CV_64F);
        this->mCostMatrix = cv::Mat::zeros(this->mCostMatrixOriginal.rows, this->mCostMatrixOriginal.cols, CV_64F);
        
        for(int row = 0; row < this->mCostMatrixOriginal.rows/*this->mRows*/; ++row) 
        {
            for(int col = 0; col < this->mCostMatrixOriginal.cols/*this->mCols*/; ++col) 
            {
                this->mCostMatrix.at<double>(row,col) = inputMatrix.at<double>(row,col);//(row < inputMatrix.rows && col < inputMatrix.cols) ? (inputMatrix.at<double>(row,col)) : std::numeric_limits<double>::max()/*0.0*/;
                this->mMaxCost = std::max(this->mMaxCost, this->mCostMatrix.at<double>(row,col));
            }
        }
        
        // Initially, none of the cells of the matrix are marked.
        this->mMarks = cv::Mat::zeros(this->mRows, this->mCols, CV_8U);
        
        if(this->mCostMatrixOriginal.rows != this->mCostMatrixOriginal.cols)
        {
            this->padMatrix(this->mMaxCost);
        }
        
        if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL) 
        {
            for(int row = 0; row < this->mRows; ++row) 
            {
                for(int col = 0; col < this->mCols; ++col) 
                {
                    this->mCostMatrix.at<double>(row,col) =  this->mMaxCost - this->mCostMatrix.at<double>(row,col);
                }
            }
        }
        else if (mode == HUNGARIAN_MODE_MINIMIZE_COST) 
        {
            // nothing to do
        }
        else 
            fprintf(stderr,"%s: unknown mode. Mode was set to HUNGARIAN_MODE_MINIMIZE_COST !\n", __FUNCTION__);
    }
    
    void Hungarian::padMatrix(double padValue)
    {
        if (padValue == 0) //KL Added check
            return;

        // we have to add some rows
        if(this->mCostMatrixOriginal.rows < this->mCostMatrixOriginal.cols)
        {
            int times = this->mCostMatrixOriginal.cols - this->mCostMatrixOriginal.rows;
            for(int i = 0; i < times; ++i)
            {
                cv::Mat row = cv::Mat::ones(1, this->mCostMatrixOriginal.cols, CV_64F);
                row *= padValue;
                this->mCostMatrix.push_back(row);
            }
        }
        else
        {
            //KL Added
            if (this->mCostMatrix.empty())
                return;

            if (this->mCostMatrix.rows < 1 )
                return;


            int times = this->mCostMatrixOriginal.rows - this->mCostMatrixOriginal.cols;
            for(int i = 0; i < times; ++i)
            {
                cv::Mat col = cv::Mat::ones(this->mCostMatrixOriginal.rows, 1, CV_64F);
                col *= padValue;
                cv::hconcat(this->mCostMatrix, col, this->mCostMatrix);
            }
        }
    }
    
    bool Hungarian::isRowCovered(int row) const
    {
        return this->mRowsCovered.at(row);
    }
    
    void Hungarian::coverRow(int row)
    {
        this->mRowsCovered.at(row) = true;
    }
    
    void Hungarian::uncoverRow(int row)
    {
        this->mRowsCovered.at(row) = false;
    }
    
    bool Hungarian::isColCovered(int col) const
    {
        return this->mColsCovered.at(col);
    }
    
    void Hungarian::coverCol(int col)
    {
        this->mColsCovered.at(col) = true;
    }
    
    void Hungarian::uncoverCol(int col)
    {
        this->mColsCovered.at(col) = false;
    }
    
    void Hungarian::clearCovers()
    {
        for(int col = 0; col < this->mCols; ++col)
        {
            this->uncoverCol(col);
        }
        
        for(int row = 0; row < this->mRows; ++row)
        {
            this->uncoverRow(row);
        }
    }
    
    bool Hungarian::isZeroStarred(int row, int col) const
    {
        return this->mMarks.at<uchar>(row, col) == STAR;
    }
    
    void Hungarian::starZero(int row, int col)
    {
        this->mMarks.at<uchar>(row, col) = STAR;
        ++this->mStarsInCol.at(col);
    }
    
    void Hungarian::unStarZero(int row, int col)
    {
        this->mMarks.at<uchar>(row, col) = NONE;
        --this->mStarsInCol.at(col);
    }
    
    bool Hungarian::colContainsStarredZero(int col) const
    {
        return this->mStarsInCol.at(col) > 0;
    }
    
    bool Hungarian::findUncoveredZero(int &zeroRow, int &zeroCol) const
    {
        for(int row = 0; row < this->mRows; ++row)
        {
            if(this->isRowCovered(row))
            {
                continue;
            }
            
            for(int col = 0; col < this->mCols; ++col)
            {
                if(this->isColCovered(col))
                {
                    continue;
                }
                
                if(this->mCostMatrix.at<double>(row, col) == 0.0)
                {
                    zeroRow = row;
                    zeroCol = col;
                    return true;
                }
            }
        }
        
        return false;
    }
    
    bool Hungarian::isZeroPrimed(int row, int col) const
    {
        return this->mMarks.at<uchar>(row, col) == PRIME;
    }
    
    void Hungarian::primeZero(int row, int col)
    {
        this->mMarks.at<uchar>(row, col) = PRIME;
    }
    
    void Hungarian::unPrimeZero(int row, int col)
    {
        this->mMarks.at<uchar>(row, col) = NONE;
    }
    
    int Hungarian::findPrimedZeroInRow(int row) const
    {
        for(int col = 0; col < this->mCols; ++col)
        {
            if(this->isZeroPrimed(row, col))
            {
                return col;
            }
        }
        
        return COLNOTFOUND;
    }
    
    void Hungarian::clearPrimes()
    {
        for(int row = 0; row < this->mRows; ++row)
        {
            for(int col = 0; col < this->mCols; ++col)
            {
                if(this->isZeroPrimed(row, col))
                {
                    this->unPrimeZero(row, col);
                }
            }
        }
    }
    
    int Hungarian::findStarInRow(int row) const
    {
        for(int col = 0; col < this->mCols; ++col)
        {
            if(this->isZeroStarred(row, col))
            {
                return col;
            }
        }
        
        return COLNOTFOUND;
    }
    
    int Hungarian::findStarInCol(int col) const
    {
        if (!this->colContainsStarredZero(col)) 
        {
            return ROWNOTFOUND;
        }
        
        for (int row = 0; row < this->mRows; ++row) 
        {
            if (this->isZeroStarred(row, col)) 
            {
                return row;
            }
        }
        
        // NOTREACHED
        return ROWNOTFOUND;
    }
    
    void Hungarian::augmentPath(const int pathLength)
    {
        int row;
        int col;
        for(int i = 0; i <= pathLength; ++i)
        {
            row = this->mPath.at(i).first;
            col = this->mPath.at(i).second;
            
            if(this->isZeroStarred(row, col))
            {
                this->unStarZero(row, col);
            }
            else
            {
                this->starZero(row, col);
            }
        }
    }
    
    double Hungarian::findSmallestUncoveredValue() const
    {
        double minValue = std::numeric_limits<double>::max();
        for(int row = 0; row < this->mRows; ++row)
        {
            if(this->isRowCovered(row))
            {
                continue;
            }
            
            for(int col = 0; col < this->mCols; ++col)
            {
                if(this->isColCovered(col))
                {
                    continue;
                }
                
                minValue = std::min(minValue, this->mCostMatrix.at<double>(row, col));
            }
        }
        
        return minValue;
    }
    
    void Hungarian::findAssignments()
    {
        this->mAssignment.clear();
        for(int row = 0; row < this->mCostMatrixOriginal.rows/*this->mRows*/; ++row)
        {
            for(int col = 0; col < this->mCostMatrixOriginal.cols/*this->mCols*/; ++col)
            {
                if(this->isZeroStarred(row, col))
                {
                    this->mAssignment.push_back(std::pair<int, int>(row, col));
                    this->mCost += this->mCostMatrixOriginal.at<double>(row, col);
                    break;
                }
            }
        }
        
        if(this->mAssignment.size() != std::min(this->mCostMatrixOriginal.rows, this->mCostMatrixOriginal.cols)/*std::min(this->mRows, this->mCols)*/)
        {
            std::cerr << std::endl << "could not provide a complete assignment" << std::endl << std::endl;
        }
    }
    
    void Hungarian::showCostMatrix()
    {
        for (int row = 0; row < this->mRows; ++row)
        {
            std::cout << "\n";
            std::cout << "     ";
            for (int col = 0; col < this->mCols; ++col)
            {
                std::cout << this->mCostMatrix.at<double>(row, col);
                if (this->isZeroStarred(row, col)) 
                {
                    std::cout << "* ";
                }
                else if (this->isZeroPrimed(row, col)) 
                {
                    std::cout << "' ";
                }
                else
                {
                    std::cout << " ";
                }
            }
        }
        std::cout << "\n";
    }
    
    void Hungarian::showMaskMatrix()
    {
        std::cout << "\n    ";
        for (int col = 0; col < this->mCols; ++col)
        {
            std::cout << " " << this->isColCovered(col);
        }
        for (int row = 0; row < this->mRows; ++row)
        {
            std::cout << "\n  " << this->isRowCovered(row) << "  ";
            for (int col = 0; col < this->mCols; ++col)
            {
                std::cout << static_cast<unsigned int>(this->mMarks.at<uchar>(row, col)) << " ";
            }
        }
        std::cout << "\n";
    }
    
    void Hungarian::showAssignment()
    {
        std::cout << "\n#(Tracker -> Point)\n[";
        for(size_t i = 0; i < this->mAssignment.size(); ++i)
        {
            if(i != this->mAssignment.size() - 1)
            {
                std::cout << "(" << this->mAssignment.at(i).first << " -> " << this->mAssignment.at(i).second << "), ";
            }
            else
            {
                std::cout << "(" << this->mAssignment.at(i).first << " -> " << this->mAssignment.at(i).second << ")";
            }
        }
        std::cout << "] => Cost = " << this->cost() << "\n";
    }
    
    /**
     * @brief Hungarian::stepOne
     * Step 1.
     * For each row of the matrix, find the smallest element and subtract it
     * from every element in its row.  Go to Step 2.
     * @return 2 => Go to Step 2.
     */
    int Hungarian::stepOne()
    {
        //KL
        if (this->mCostMatrix.empty())
            return 7; //Jump to step 7?

        for(int row = 0; row < this->mRows; ++row)
        {
            double minCostInRow = this->mCostMatrix.at<double>(row, 0);
            for(int col = 1; col < this->mCols; ++col)
            {
                minCostInRow = std::min(minCostInRow, this->mCostMatrix.at<double>(row, col));
            }
            
            for(int col = 0; col < this->mCols; ++col)
            {
                this->mCostMatrix.at<double>(row, col) -= minCostInRow;
            }
        }
        
        return 2;
    }
    
    /**
     * @brief Hungarian::stepTwo
     * Find a zero (Z) in the resulting matrix.  If there is no starred zero 
     * in its row or column, star Z. Repeat for each element in the matrix. 
     * Go to Step 3.
     * @return 
     */
    int Hungarian::stepTwo()
    {
        //KL
        if (this->mCostMatrix.empty())
            return 7; //Jump to step 7?

        for(int row = 0; row < this->mRows; ++row)
        {
            if(this->isRowCovered(row))
            {
                continue;
            }
            
            for(int col = 0; col < this->mCols; ++col)
            {
                if(this->isColCovered(col))
                {
                    continue;
                }
                
                if(this->mCostMatrix.at<double>(row,col) == 0.0)
                {
                    this->starZero(row,col);
                    this->coverRow(row);
                    this->coverCol(col);
                    break;
                }
            }
        }
        this->clearCovers();
        return 3;
    }
    
    /**
     * @brief Hungarian::stepThree
     * Cover each column containing a starred zero. If all columns are covered, 
     * the starred zeros describe a complete set of unique assignments.  
     * In this case, terminate the algorithm, otherwise, Go to Step 4.
     * @return 
     */
    int Hungarian::stepThree()
    {
        int coveredColums = 0;
        
        for(int row = 0; row < this->mRows; ++row)
        {
            for(int col = 0; col < this->mCols; ++col)
            {
                if(this->isZeroStarred(row, col))
                {
                    this->coverCol(col);
                }
            }
        }
        
        for(int col = 0; col < this->mCols; ++col)
        {
            if(this->isColCovered(col))
            {
                ++coveredColums;
            }
        }
        
        if(coveredColums >= this->mCols || coveredColums >= this->mRows)
        {
            return 7;
        }
        
        return 4;
    }
    
    /**
     * @brief Hungarian::stepFour
     * Find a noncovered zero and prime it. If there is no starred zero in 
     * the row containing this primed zero, Go to Step 5. Otherwise, cover 
     * this row and uncover the column containing the starred zero. 
     * Continue in this manner until there are no uncovered zeros left. 
     * Save the smallest uncovered value and Go to Step 6.
     * @return 
     */
    int Hungarian::stepFour()
    {
        int zeroRow = -1;
        int zeroCol = -1;
        int starCol = COLNOTFOUND;
        while(true)
        {
            if(!this->findUncoveredZero(zeroRow, zeroCol))
            {
                return 6;
            }
            else
            {
                this->primeZero(zeroRow, zeroCol);
                starCol = this->findStarInRow(zeroRow);
                if(starCol != COLNOTFOUND)
                {
                    this->coverRow(zeroRow);
                    this->uncoverCol(starCol);
                }
                else
                {
                    this->mPath.at(0).first = zeroRow;
                    this->mPath.at(0).second = zeroCol;
                    return 5;
                }
            }
        }
    }
    
    /**
     * @brief Hungarian::stepFive
     * Construct a series of alternating primed and starred zeros as follows. 
     * Let Z0 represent the uncovered primed zero found in Step 4. Let Z1 denote 
     * the starred zero in the column of Z0 (if any). Let Z2 denote the primed 
     * zero in the row of Z1 (there will always be one). Continue until the 
     * series terminates at a primed zero that has no starred zero in its column. 
     * Unstar each starred zero of the series, star each primed zero of the series, 
     * erase all primes and uncover every line in the matrix. 
     * Return to Step 3.
     * @return 
     */
    int Hungarian::stepFive()
    {
        bool done       = false;
        int row         = -1;
        int col         = -1;
        int pathLength   = 0;
        
        while(!done)
        {
            // First construct the alternating path...
            row = findStarInCol(this->mPath.at(pathLength).second);
            
            if(row != ROWNOTFOUND)
            {
                ++pathLength;
                this->mPath.at(pathLength).first = row;
                this->mPath.at(pathLength).second = this->mPath.at(pathLength-1).second;
            }
            else
            {
                done = true;
            }
            
            if(!done)
            {
                col = this->findPrimedZeroInRow(this->mPath.at(pathLength).first);
                ++pathLength;
                this->mPath.at(pathLength).first = this->mPath.at(pathLength-1).first;
                this->mPath.at(pathLength).second = col;
            }
        }
        
        // Then modify it.
        this->augmentPath(pathLength);
        
        this->clearCovers();
        this->clearPrimes();
        
        return 3;
    }
    
    /**
     * @brief Hungarian::stepSix
     * Add the smallest uncovered value in the costmatrix to every element of each 
     * covered row, and subtract it from every element of each uncovered column.
     * Return to Step 4 without altering any stars, primes, or covered lines.
     * @return 
     */
    int Hungarian::stepSix()
    {
        double minValue = this->findSmallestUncoveredValue();
        
        for(int row = 0; row < this->mRows; ++row)
        {
            for(int col = 0; col < this->mCols; ++col)
            {
                if(this->isRowCovered(row))
                {
                    this->mCostMatrix.at<double>(row, col) += minValue;
                }
                
                if(!this->isColCovered(col))
                {
                    this->mCostMatrix.at<double>(row, col) -= minValue;
                }
            }
        }
        
        return 4;
    }
    
    /**
     * @brief Hungarian::stepSeven
     */
    void Hungarian::stepSeven()
    {
        this->findAssignments();
        if(DEBUGMODE)
        {
            std::cout << "---------Run Complete----------";
            this->showCostMatrix();
            this->showMaskMatrix();
            this->showAssignment();
            std::cout << "-------------------------------" << std::endl;
        }
//        else
//        {
//            std::cout << "---------Run Complete----------";
//            this->showAssignment();
//            std::cout << "-------------------------------" << std::endl;
//        }
    }
    
    void Hungarian::runMunkres()
    {
        if(this->mIsInit)
        {
            bool done = false;
            int step = 1;
            while(!done)
            {
                if(DEBUGMODE)
                {
                    std::cout << "\n";
                    std::cout << "------------Step {" << step << "}-------------";
                    std::cout <<"\n";
                    this->showCostMatrix();
                    this->showMaskMatrix();
                    std::cout <<"\n";
                }
                
                switch (step) 
                {
                    case 1:
                        step = this->stepOne();
                        break;
                    case 2:
                        step = this->stepTwo();
                        break;
                    case 3:
                        step = this->stepThree();
                        break;
                    case 4:
                        step = this->stepFour();
                        break;
                    case 5:
                        step = this->stepFive();
                        break;
                    case 6:
                        step = this->stepSix();
                        break;
                    case 7:
                        this->stepSeven();
                        done = true;
                        break;
                }
            }
        }
    }
    
    std::vector<std::pair<int, int> > const& Hungarian::getAssignment() const
    {
        return this->mAssignment;
    }
    
    cv::Mat Hungarian::getAssignmentAsMatrix()
    {
        cv::Mat ret = this->mMarks(cv::Rect(0, 0 , this->mCostMatrixOriginal.cols, this->mCostMatrixOriginal.rows)).clone();
        return ret;
    }
    
}

