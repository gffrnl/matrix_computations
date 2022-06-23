#ifndef MATRIX_COMPUTATIONS_MATRIX_H_INCLUDED
#define MATRIX_COMPUTATIONS_MATRIX_H_INCLUDED

#include <iostream>
#include <valarray>
#include <algorithm>
#include <cmath>
#include <limits>


/*
 * helpers: output operator << for std::valarray, std::slice_array and std::gslice_array
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& va)
{
    auto it = std::begin(va);
    for (; it != std::end(va)-1; ++it)
        os << *it << "  ";
    os << *it;
    return os;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::slice_array<T>& sa)
{
    return os << std::valarray<T>(sa);
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::gslice_array<T>& ga)
{
    return os << std::valarray<T>(ga);
}


/*
 * In the following we define template classes to representing general matrices
 *
 */

namespace matrix_computations {

    template<typename T>
    class Matrix { // FORTRAN-like storage and indexing
    // References:
    // https://en.cppreference.com/w/cpp/numeric/valarray/slice
    // https://en.cppreference.com/w/cpp/numeric/valarray/slice_array
    protected:
        std::valarray<T> data;
        std::size_t dim;

    public:
        Matrix(std::size_t nrows, std::size_t ncols, T init = 0)
            : data(init, nrows*ncols), dim(nrows) {}
        virtual ~Matrix() {}

        std::size_t numb_entries() const { return data.size(); }
        std::size_t numb_rows()    const { return dim; }
        std::size_t numb_columns() const { return data.size()/numb_rows(); }

        bool is_square() const { return numb_rows() == numb_columns(); }
        bool is_null() const
        {
            if (std::abs(data).max() > std::numeric_limits<T>::epsilon())
                return false;
            else
                return true;
        }

        T& operator()(std::size_t i, std::size_t j)
        {
            return data[(i-1) + (j-1)*numb_rows()];
        }

        std::ostream& print(std::ostream& os) const // printing row by row
        {
            os << "{ ";
            for (std::size_t i = 1; i < numb_rows(); ++i) {
                os << "{";
                std::valarray<T> r = row(i);
                auto it = std::begin(r);
                for (; it != std::end(r)-1; ++it)
                    os << *it << ",  ";
                os << *it;
                os << "}, ";
            }
            os << "{";
            std::valarray<T> r = row(numb_rows());
            auto it = std::begin(r);
            for (; it != std::end(r)-1; ++it)
                os << *it << ",  ";
            os << *it;
            os << "}";
            os << " }";
            return os;
        }


    //    std::slice_array<T> row(std::size_t i)
    //    {
    //        return data[std::slice(i-1, numb_columns(), numb_rows())];
    //    }
    //    std::slice_array<T> column(std::size_t j)
    //    {
    //        return data[std::slice((j-1)*numb_rows(), numb_rows(), 1)];
    //    }

        std::valarray<T> row(std::size_t i) const
        {
            return data[std::slice(i-1, numb_columns(), numb_rows())];
        }
        std::valarray<T> column(std::size_t j) const
        {
            return data[std::slice((j-1)*numb_rows(), numb_rows(), 1)];
        }

        void swap_rows(std::size_t i1, std::size_t i2)
        {
            const std::valarray<T> temp = row(i1);
            row(i1) = row(i2);
            row(i2) = temp;
        }
        void swap_columns(std::size_t j1, std::size_t j2)
        {
            const std::valarray<T> temp = column(j1);
            column(j1) = column(j2);
            column(j2) = temp;
        }

        void clear(T value = 0) { data = value; }

        void transpose()
        {
            std::valarray<T> temp =
                data[std::gslice(0, {numb_rows(), numb_columns()}, {1, numb_rows()})];
            dim = numb_columns();
            data = temp;
        }

        void augment(const Matrix<T>& m)
        {
            if (numb_rows() != m.numb_rows())
                throw std::invalid_argument("the original matrix and the passed matrix "
                                            "have different number of rows");
            std::valarray<T> temp = data;
            data.resize(temp.size() + m.data.size());
            auto it  = std::begin(data);
            auto itt = std::begin(temp);
            auto itm = std::begin(m.data);
            for(; itt != std::end(temp); ++it, ++itt)
                *it = *itt;
            for(; itm != std::end(m.data); ++it, ++itm)
                *it = *itm;
        }



        /*
         *  Gaussian eliminations
         */
        void gauss_elimination_piv_none()
        {
            T pivot, ratio;
            std::size_t h = 1; // initialization of the pivot row
            std::size_t k = 1; // initialization of the pivot column
            while (h < numb_rows()  &&  k < numb_columns()) {
                pivot = (*this)(h, k);
                // do for all rows below pivot
                for (std::size_t i = h+1; i <= numb_rows(); ++i) {
                    ratio = (*this)(i, k) / pivot;
                    // fill with zeros the lower part of pivot column
                    (*this)(i, k) = 0;
                    for (std::size_t j = k+1; j <= numb_columns(); ++j)
                        (*this)(i, j) -= (*this)(h, j) * ratio;
                }
                ++h;
                ++k;
            }
        }
        void gauss_elimination_piv_partial()
        {
            T max, pivot, ratio;
            std::size_t i_max;
            std::size_t h = 1; // initialization of the pivot row
            std::size_t k = 1; // initialization of the pivot column
            while (h < numb_rows()  &&  k < numb_columns()) {
                // find the k-th pivot
                max = 0;
                i_max = 0;
                for (std::size_t i = h; i <= numb_rows(); ++i)
                    if (std::abs((*this)(i, k)) > max)
                        i_max = i;

                if (i_max == 0) // no pivot in this column, pass to next column
                        ++k;
                else {
                    swap_rows(i_max, h);
                    pivot = (*this)(h, k);
                    // do for all rows below pivot
                    for (std::size_t i = h+1; i <= numb_rows(); ++i) {
                        ratio = (*this)(i, k) / pivot;
                        // fill with zeros the lower part of pivot column
                        (*this)(i, k) = 0;
                        for (std::size_t j = k+1; j <= numb_columns(); ++j)
                            (*this)(i, j) -= (*this)(h, j) * ratio;
                    }
                    // increase pivot row and column
                    ++h;
                    ++k;
                }
            }
        }
        void gauss_elimination_piv_complete() // REVISAR!!
        {
            if (is_null()) // null matrix, return to the caller
                return;
            T max, pivot, ratio;
            std::size_t i_max, j_max;
            std::size_t h = 1; // initialization of the pivot row
            std::size_t k = 1; // initialization of the pivot column
            while (h < numb_rows()  &&  k < numb_columns()) {
                // find the pivot
                max = 0;
                i_max = 0;
                j_max = 0;
                for (std::size_t i = h; i <= numb_rows(); ++i)
                    for (std::size_t j = k; j < numb_columns(); ++j)
                        if (std::abs((*this)(i, k)) > max) {
                            i_max = i;
                            j_max = j;
                        }

                swap_columns(j_max, k);
                swap_rows(i_max, h);
                pivot = (*this)(h, k);
                // do for all rows below pivot
                for (std::size_t i = h+1; i <= numb_rows(); ++i) {
                    ratio = (*this)(i, k) / pivot;
                    // fill with zeros the lower part of pivot column
                    (*this)(i, k) = 0;
                    for (std::size_t j = k+1; j <= numb_columns(); ++j)
                        (*this)(i, j) -= (*this)(h, j) * ratio;
                }
                // increase pivot row and column
                ++h;
                ++k;
            }
        }
        void gauss_elimination()
        {
            this->gauss_elimination_piv_partial();
        }

    }; // end template class template<typename T> matrix_computations::Matrix<T>



    template<typename T>
    class SquareMatrix : public Matrix<T> {
    public:
        SquareMatrix(std::size_t order, T init = 0)
            : Matrix<T>(order, order, init) {}
        virtual ~SquareMatrix() {}

    private:
	using Matrix<T>::augment;

    public:
        std::size_t order() const { return Matrix<T>::numb_rows(); }

    //    void swap_row_and_column(std::size_t i, std::size_t j)
    //    {
    //        const std::valarray<T> temp = Matrix<T>::row(i);
    //        Matrix<T>::row(i) = Matrix<T>::column(j);
    //        Matrix<T>::column(j) = temp;
    //    }

        void identify(T value = 1)
        {
            Matrix<T>::clear();
            diagonal() = value;
        }

    //    std::slice_array<T> diagonal()
    //    {
    //        return Matrix<T>::data[std::slice(0, order(), order()+1)];
    //    }
    //    std::slice_array<T> secondary_diagonal()
    //    {
    //        return Matrix<T>::data[std::slice(order()-1, order(), order()-1)];
    //    }

        std::valarray<T> diagonal() const
        {
            return Matrix<T>::data[std::slice(0, order(), order()+1)];
        }
        std::valarray<T> secondary_diagonal() const
        {
            std::valarray<T> temp = Matrix<T>::data[std::slice(order()-1, order(), order()-1)];
            std::reverse(std::begin(temp), std::end(temp));
            return temp;
        }


        T trace() const
        {
            return std::valarray<T>(diagonal()).sum();
        }
        T determinant() const
        {
            // TODO
            throw std::runtime_error("`determinant()` not implemented");
            return 0;
        }
    }; // end template class template<typename T> matrix_computations::SquareMatrix<T>



    namespace solvers {
        template<typename T>
        void gauss_jordan(const Matrix<T>& A, const Matrix<T>& b)
        {

        }
    } // end namespace matrix_computations::solvers


} // end namespace matrix_computations


/*
 * non-member methods
 *
 */

template<typename T>
std::ostream& operator<<(std::ostream& os, const matrix_computations::Matrix<T>& m)
{
    return m.print(os);
}

#endif // MATRIX_COMPUTATIONS_MATRIX_H_INCLUDED
