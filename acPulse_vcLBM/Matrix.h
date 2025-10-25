#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
	public:
			// Constructor
			Matrix() {matrix_=NULL;}

			Matrix(unsigned short lines, unsigned short columns) //lines: no of lines; columns: no. of columns
			{
				lines_=lines;
				columns_=columns;
				matrix_= new double*[lines_];
				for (unsigned short i=0;i<lines_;i++) matrix_[i]=new double[columns_];
				//initialize elements
				for (unsigned short i=0;i<lines_;i++)
					for(unsigned short j=0;j<columns_;j++)
						matrix_[i][j]=0.;
			}

			// Copy Constructor
			Matrix(const Matrix &matrix)
			{
				lines_=matrix.getNoOfLines();
				columns_=matrix.getNoOfColumns();
				matrix_= new double*[lines_];
				for (unsigned short i=0;i<lines_;i++) matrix_[i]=new double[columns_];

				for (unsigned short i=0;i<lines_;i++)
					for(unsigned short j=0;j<columns_;j++)
						matrix_[i][j]=matrix.getElement(i,j);
			}

			// Assignment operator
			Matrix& operator= (const Matrix& matrix)
			{
				lines_=matrix.getNoOfLines();
				columns_=matrix.getNoOfColumns();
				for (unsigned short i=0;i<lines_;i++)
					for(unsigned short j=0;j<columns_;j++)
						matrix_[i][j]=matrix.getElement(i,j);

				return *this;
			}

			Matrix(Matrix* matrix)
			{
				lines_=matrix->getNoOfLines();
				columns_=matrix->getNoOfColumns();
				matrix_= new double*[lines_];
				for (unsigned short i=0;i<lines_;i++) matrix_[i]=new double[columns_];

				for (unsigned short i=0;i<lines_;i++)
					for(unsigned short j=0;j<columns_;j++)
						matrix_[i][j]=matrix->getElement(i,j);
			}

			
			//destructor
			~Matrix()
			{
				if (matrix_)
				{
					for (unsigned short i=0;i<lines_;i++) delete[] matrix_[i];
					delete[] matrix_;
				}
			}

	// Members
			double getElement(unsigned short line,unsigned short column) const 
			{
				if ((line<lines_) && (column<columns_)) return matrix_[line][column];
				std::cout << "\nMatrix element in line " << line << " column " << column << " does not exist. This is a" << lines_ << "x" << columns_ << " matrix."; 
				return 0.;
			}
			Matrix* getMatrix() {return this;}

			unsigned short getNoOfLines() const  {return lines_;}
			unsigned short getNoOfColumns() const {return columns_;}

			void setElement(unsigned short line,unsigned short column,double value)
			{	
				if ((line<lines_) && (column<columns_)) matrix_[line][column]=value;
				else std::cout << "\nMatrix element in line " << line << " column " << column << " does not exist.This is a" << lines_ << "x" << columns_ << " matrix.";
			}


	private:
			unsigned short columns_;
			unsigned short lines_;
			double** matrix_;
};

#endif