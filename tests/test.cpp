/* Library matrix_computations unit tests */

#include "matrix_computations/Matrix.h"
#include <iostream>

using namespace matrix_computations;
using namespace std;

template<typename T>
void info(std::string label, const Matrix<T>& m)
{
    cout << " Matrix " << label << '\n' << endl;
    cout << "\t" << label << " = " << m << '\n' << endl;
    cout << "\tProperties:" << endl;
    cout << "\tnumber of entries: " << m.numb_entries() << endl;
    cout << "\tnumber of rows:    " << m.numb_rows()    << endl;
    cout << "\tnumber of columns: " << m.numb_columns() << endl;
    cout << "\tis_square?:        " << boolalpha << m.is_square()    << endl;
    cout << endl;
    for (std::size_t i = 1; i <= m.numb_rows(); ++i)
        cout << "\trow    " << i << ": " << m.row(i) << endl;
    cout << endl;
    for (std::size_t j = 1; j <= m.numb_columns(); ++j)
        cout << "\tcolumn " << j << ": " << m.column(j) << endl;
}

template<typename T>
void test_swap_columns(std::string label, Matrix<T> m, std::size_t i1, std::size_t i2)
{
    cout << "Swapping columns " << i1 << " and " << i2 << " of matrix " << label << "..." << endl;
    m.swap_columns(i1, i2);
    info(label, m);
    cout << '\n' << endl;
}

template<typename T>
void test_swap_rows(std::string label, Matrix<T> m, std::size_t i1, std::size_t i2)
{
    cout << "Swapping rows " << i1 << " and " << i2 << " of matrix " << label << "..." << endl;
    m.swap_rows(i1, i2);
    info(label, m);
    cout << '\n' << endl;
}

/*
template<typename T>
void test_swap_row_and_column(std::string label, SquareMatrix<T> m, std::size_t i, std::size_t j)
{
    cout << "Swapping row" << i << " and column" << j << " of square matrix " << label << "..." << endl;
    m.swap_row_and_column(i, j);
    info(label, m);
    cout << '\n' << endl;
}
*/

template<typename T>
void test_transpose(std::string label, Matrix<T> m)
{
    cout << "Transposing  matrix " << label << "..." << endl;
    m.transpose();
    info(label, m);
    cout << '\n' << endl;
}

template<typename T>
void test_identify(std::string label, SquareMatrix<T> m, T val = 1)
{
    cout << "Identifying square matrix " << label << " to diagonal value " << val << "..." << endl;
    m.identify(val);
    info(label, m);
    cout << '\n' << endl;
}

template<typename T>
void test_gauss_elimination_piv_none(std::string label, Matrix<T> m)
{
    cout << "Performing Gauss elimination without pivoting on matrix " << label << "..." << endl;
    m.gauss_elimination_piv_none();
    info(label, m);
    cout << '\n' << endl;
}
template<typename T>
void test_gauss_elimination_piv_partial(std::string label, Matrix<T> m)
{
    cout << "Performing Gauss elimination with partial pivoting on matrix " << label << "..." << endl;
    m.gauss_elimination();
    info(label, m);
    cout << '\n' << endl;
}
template<typename T>
void test_gauss_elimination_piv_complete(std::string label, Matrix<T> m)
{
    cout << "Performing Gauss elimination with complete pivoting on matrix " << label << "..." << endl;
    m.gauss_elimination_piv_complete();
    info(label, m);
    cout << '\n' << endl;
}

template<typename T>
void test_augment(std::string label, Matrix<T> m, std::string label2, const Matrix<T>& m2)
{
    cout << "Augmenting matrix " << label << " with matrix " << label2 << "..." << endl;
    m.augment(m2);
    info(label, m);
    cout << '\n' << endl;
}


int main(int argc, char* argv[])
{
    // O que testar?
    //
    SquareMatrix<double> m1(3);
    m1(1,1)=1; m1(1,2)= 3; m1(1,3)= 1;
    m1(2,1)=1; m1(2,2)= 1; m1(2,3)=-1;
    m1(3,1)=3; m1(3,2)=11; m1(3,3)= 5;

    Matrix<double> m1_4(3, 1);
    m1_4(1,1)= 9;
    m1_4(2,1)= 1;
    m1_4(3,1)=35;

    cout << '\n';
    info("m1", m1);
    cout << '\n' << endl;
    info("m1_4", m1_4);
    cout << '\n' << endl;

    //test_swap_columns("m1", m1, 1, 2);
    //test_swap_columns("m1", m1, 1, 3);
    //test_swap_columns("m1", m1, 2, 3);

    //test_swap_rows("m1", m1, 1, 2);
    //test_swap_rows("m1", m1, 1, 3);
    //test_swap_rows("m1", m1, 2, 3);

    //test_swap_rows("m1_4", m1_4, 1, 2);
    //test_swap_rows("m1_4", m1_4, 1, 3);
    //test_swap_rows("m1_4", m1_4, 2, 3);

////    //test_swap_row_and_column("m1", m1, 1, 1);
////    //test_swap_row_and_column("m1", m1, 1, 2);
////    //test_swap_row_and_column("m1", m1, 1, 3);

    //test_identify("m1", m1, 15.44);

    //test_transpose("m1", m1);
    //test_transpose("m1_4", m1_4);

    //cout << "order of m1:              " << m1.order()              << endl;
    //cout << "diagonal of m1:           " << m1.diagonal()           << endl;
    //cout << "secondary diagonal of m1: " << m1.secondary_diagonal() << endl;
    //cout << "trace of m1:              " << m1.trace()              << endl;
    //cout << "determinant of m1:        " << m1.determinant()        << endl;

    //test_gauss_elimination_piv_none("m1", m1);
    //test_gauss_elimination_piv_partial("m1", m1);
    //test_gauss_elimination_piv_complete("m1", m1);

    //test_augment("m1", m1, "m1_4", m1_4);

    Matrix<double> m1_aug = m1;
    m1_aug.augment(m1_4);

    test_gauss_elimination_piv_none("m1_aug", m1_aug);
    test_gauss_elimination_piv_partial("m1_aug", m1_aug);
    test_gauss_elimination_piv_complete("m1_aug", m1_aug);

    return 0;
}
