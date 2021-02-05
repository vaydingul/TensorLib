#include <random>
#include <chrono>
#include <iostream>

template <typename T>
Tensor<T>::Tensor()
{
    /*
It is the default constructor. It creates one dim1 and dim2umn ( a point actually) ; then, it sets to zero.
*/
}

template <typename T>
Tensor<T>::Tensor(const int dim1, const int dim2, const int dim3, const std::string TYPE)
{
    /*
This constructor takes necessary dim1 and dim2umn number as input and creates 2d std::vector with this information.
Then, it sets all the elements of this std::vector according to the third input argument. It decides which type of information
the element of this Tensor should convey. The TYPE argument can be "NORMAL" or "ZEROS".

Example:
Tensor<double> A(5,5,"NORMAL");
*/
    this->dim1 = dim1;
    this->dim2 = dim2;
    this->dim3 = dim3;
    this->self.resize(this->dim1);
    for (int i = 0; i < this->dim1; i++)
    {
        this->self[i].resize(this->dim2);

        for (int j = 0; j < this->dim3; j++)
        {
            this->self[i][j].resize(this->dim3);
        }
    }

    if (TYPE == "NORMAL")
    {
        int seed = chrono::system_clock::now().time_since_epoch().count();
        mt19937 engine(seed);
        normal_distribution<T> distribution(0.0, 1.0);

        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = distribution(engine);
                }
            }
        }
    }
    else if (TYPE == "ZEROS")
    {
        
        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = 0.0;
                }
            }
        }
    }

    else if (TYPE == "ONES")
    {
        
        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = 1.0;
                }
            }
        }
    }
}


template <typename T>
Tensor<T>::Tensor(const int dim1, const int dim2, const std::string TYPE)
// TODO: Add more option. E.g. "ONES" for identity Tensor.
{
    /*
This constructor takes necessary dim1 and dim2umn number as input and creates 2d std::vector with this information.
Then, it sets all the elements of this std::vector according to the third input argument. It decides which type of information
the element of this Tensor should convey. The TYPE argument can be "NORMAL" or "ZEROS".

Example:
Tensor<double> A(5,5,"NORMAL");
*/
    this->dim1 = dim1;
    this->dim2 = dim2;
    this->dim3 = 1;
    this->self.resize(this->dim1);
    for (int i = 0; i < this->dim1; i++)
    {
        this->self[i].resize(this->dim2);

        for (int j = 0; j < this->dim3; j++)
        {
            this->self[i][j].resize(this->dim3);
        }
    }

    if (TYPE == "NORMAL")
    {
        int seed = chrono::system_clock::now().time_since_epoch().count();
        mt19937 engine(seed);
        normal_distribution<T> distribution(0.0, 1.0);

        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = distribution(engine);
                }
            }
        }
    }
    else if (TYPE == "ZEROS")
    {
        
        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = 0.0;
                }
            }
        }
    }

    else if (TYPE == "ONES")
    {
        
        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = 1.0;
                }
            }
        }
    }
}

template <typename T>
Tensor<T>::Tensor(const int dim1, const std::string TYPE)
// TODO: Add more option. E.g. "ONES" for identity Tensor.
{
    /*
This constructor takes necessary dim1 and dim2umn number as input and creates 2d std::vector with this information.
Then, it sets all the elements of this std::vector according to the third input argument. It decides which type of information
the element of this Tensor should convey. The TYPE argument can be "NORMAL" or "ZEROS".

Example:
Tensor<double> A(5,5,"NORMAL");
*/
    this->dim1 = dim1;
    this->dim2 = 1;
    this->dim3 = 1;
    this->self.resize(this->dim1);
    for (int i = 0; i < this->dim1; i++)
    {
        this->self[i].resize(this->dim2);

        for (int j = 0; j < this->dim3; j++)
        {
            this->self[i][j].resize(this->dim3);
        }
    }

    if (TYPE == "NORMAL")
    {
        int seed = chrono::system_clock::now().time_since_epoch().count();
        mt19937 engine(seed);
        normal_distribution<T> distribution(0.0, 1.0);

        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = distribution(engine);
                }
            }
        }
    }
    else if (TYPE == "ZEROS")
    {
        
        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = 0.0;
                }
            }
        }
    }

    else if (TYPE == "ONES")
    {
        
        for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = 1.0;
                }
            }
        }
    }
}


template <typename T>
Tensor<T>::Tensor(const Tensor<T> &rhs)
{
    /*
Default copy constructor.

Example:
Tensor<double> A(B);
*/
    this->dim1 = rhs.getSize()(0);
    this->dim2 = rhs.getSize()(1);
    this->dim3 = rhs.getSize()(2);

    this->self.resize(this->dim1);
    for (int i = 0; i < this->dim1; i++)
    {
        this->self[i].resize(this->dim2);

        for (int j = 0; j < this->dim3; j++)
        {
            this->self[i][j].resize(this->dim3);
        }
    }

     for (int i = 0; i < this->dim1; i++)
        {
            for (int j = 0; j < this->dim2; j++)
            {
                for (int k = 0; k < this->dim3; k++)
                {
                    this->self[i][j][k] = rhs.self[i][j][k];
                }
            }
        }
}

template <typename T>
Tensor<T>::~Tensor()
{
    /*
Default deconstructor.
*/
}
template <typename T>
Tensor<T> Tensor<T>::concat(const Tensor<T> &A, const int direction)
{   // ! MAKE MORE ROBUST.

    
    if (direction == "ver")
    {
        Tensor<T> concatTensor(this->dim1 + A.getSize()(0), this->dim2, "ZEROS");

        for (int m = 0; m < this->dim2; m++)
        {
            for (int i = 0; i < this->dim1; i++)
            {
                concatTensor(i, m) = this->self[i][m];
            }
            for (int i = 0; i < A.getSize()(0); i++)
            {
                concatTensor(i + this->dim1, m) = A(i, m);
            }
        }
        return concatTensor;
    }
    else if (direction == "hor")
    {
        Tensor<T> concatTensor(this->dim1, this->dim2 + A.getSize()(1), "ZEROS");

        for (int m = 0; m < this->dim1; m++)
        {
            for (int i = 0; i < this->dim2; i++)
            {
                concatTensor(m, i) = this->self[m][i];
            }
            for (int i = 0; i < A.getSize()(1); i++)
            {
                concatTensor(m, i + this->dim2) = A(m, i);
            }
        }
        return concatTensor;
    }
}

template <typename T>
Tensor<T> Tensor<T>::operator=(const Tensor<T> &rhs)
{
    /*
Copy assignment operator.

Example:
A = B;
*/

    this->dim1 = rhs.dim1;
    this->dim2 = rhs.dim2;

    self.resize(this->dim1);
    for (int i = 0; i < this->self.size(); i++)
    {
        this->self[i].resize(this->dim2);
    }

    for (int i = 0; i < this->dim1; i++)
    {
        for (int j = 0; j < this->dim2; j++)
        {
            this->self[i][j] = rhs.self[i][j];
        }
    }

    return *this;
}

/* Tensor Tensor<T>::operator=(Tensor &rhs)
{

    int new_dim1s = rhs.getSize()(0);
    int new_dim2s = rhs.getSize()(1);

    self.resize(new_dim1s);
    for (int i = 0; i < self.size(); i++)
    {
        this->self[i].resize(new_dim2s);
    }

    for (int i = 0; i < new_dim1s; i++)
    {
        for (int j = 0; j < new_dim2s; j++)
        {
            this->self[i][j] = rhs.self[i][j];
        }
    } 

    this->dim1 = new_dim1s;
    this->dim2 = new_dim2s;

    return *this;
}*/
/*
// Addition of Two Matrices
Tensor Tensor<T>::operator+(Tensor &B)
{
    Tensor sum(this->dim1, this->dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            sum(i, j) = this->self[i][j] + B(i, j);
        }
    }
    return sum;
}*/
template <typename T>
Tensor<T> Tensor<T>::operator+(const Tensor<T> &B)
{
    /* 
    Addition operator for Tensor class. It requires same shape of matrices to be able to manage the process.

    Example:
    C = A + B;
    */
    Tensor sum(this->dim1, this->dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            sum(i, j) = this->self[i][j] + B(i, j);
        }
    }
    return sum;
}

/* Tensor Tensor<T>::operator+=(Tensor &B)
{
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            this->self[i][j] += B(i, j);
        }
    }
    return *this;
} */
template <typename T>
Tensor<T> Tensor<T>::operator+=(const Tensor<T> &B)
{
    /* 
    Addition and Assign operator for Tensor class. It requires same shape of matrices to be able to manage the process.

    Example:
    C += A;
    */
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            this->self[i][j] += B(i, j);
        }
    }
    return *this;
}

/* // Subtraction of Two Matrices
Tensor Tensor<T>::operator-(Tensor &B)
{
    Tensor diff(this->dim1, this->dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            diff(i, j) = this->self[i][j] - B(i, j);
        }
    }

    return diff;
} */

// Subtraction of Two Matrices
template <typename T>
Tensor<T> Tensor<T>::operator-(const Tensor<T> &B)
{
    /* 
    Substraction operator for Tensor class. It requires same shape of matrices to be able to manage the process.

    Example:
    C = A - B;
    */
    Tensor diff(this->dim1, this->dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {

            diff(i, j) = this->self[i][j] - B(i, j);
        }
    }

    return diff;
}

/* Tensor Tensor<T>::operator-=(Tensor &B)
{
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            this->self[i][j] -= B(i, j);
        }
    }

    return *this;
} */

// Subtraction of Two Matrices
template <typename T>
Tensor<T> Tensor<T>::operator-=(const Tensor<T> &B)
{
    /* 
    Substraction and Assign operator for Tensor class. It requires same shape of matrices to be able to manage the process.

    Example:
    C -= A ;
    */
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {

            this->self[i][j] -= B(i, j);
        }
    }

    return *this;
}

/* Tensor Tensor<T>::operator*(Tensor &B)
{
    Tensor multip(dim1, B.getSize()(1), 0.0);

    int i, j, k;
    for (i = 0; i < this->dim1; i++)
    {
        for (k = 0; k < this->dim2; k++)
        {
            for (j = 0; j < B.getSize()(1); j++)
            {
                multip(i, j) += this->self[i][k] * B(k, j);
            }
            //std::cout << multip(i,j) << " ";
        }
        //std::cout << endl;
    }
    return multip;
} */
template <typename T>
Tensor<T> Tensor<T>::operator*(const Tensor<T> &B)
{
    /* 
    Multiplication operator for Tensor class. It requires the fact that the dim2umn number of (lHS) Tensor should be equal to the 
    dim1 number of (RHS) Tensor. 

    Example:
    C = A * B;
    */
    Tensor multip(dim1, B.getSize()(1), 0.0);

    int i, j, k;
    for (i = 0; i < this->dim1; i++)
    {
        for (k = 0; k < this->dim2; k++)
        {
            for (j = 0; j < B.getSize()(1); j++)
            {
                multip(i, j) += this->self[i][k] * B(k, j);
            }
            //std::cout << multip(i,j) << " ";
        }
        //std::cout << endl;
    }
    return multip;
}

/* Tensor Tensor<T>::operator*=(Tensor &B)
{
    Tensor multip;

    multip = (*this) * B;
    (*this) = multip;
    return *this;
} */
template <typename T>
Tensor<T> Tensor<T>::operator*=(const Tensor<T> &B)
{
    /* 
    Left Multiplication and Assign operator for Tensor class.  

    Example:
    C *= A;

    Above example actually does the following:

    Tensor<double> C,A;

    C *= A;   =====>  C = C * A; 
    */
    Tensor multip;

    multip = (*this) * B;
    (*this) = multip;
    return *this;
}
template <typename T>
Tensor<T> Tensor<T>::operator+(T scalar)
{
    /* 
    It is a scalar addition operator. It performs an addition operation on every element in the Tensor ( with same value).

    Example:

    A = B + 2;
    */

    Tensor result(this->dim1, this->dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            result(i, j) = this->self[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
template <typename T>
Tensor<T> Tensor<T>::operator-(T scalar)
{
    /* 
    It is a scalar substraction operator. It performs an substraction operation on every element in the Tensor ( with same value).

    Example:

    A = B - 2;
    */
    Tensor result(dim1, dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            result(i, j) = this->self[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
template <typename T>
Tensor<T> Tensor<T>::operator*(T scalar)
{
    /* 
    It is a scalar multiplication operator. It performs an multiplication operation on every element in the Tensor ( with same value).

    Example:

    A = B * 2;
    */
    Tensor result(dim1, dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            result(i, j) = this->self[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
template <typename T>
Tensor<T> Tensor<T>::operator/(T scalar)
{
    /* 
    It is a scalar division operator. It performs an division operation on every element in the Tensor ( with same value).

    Example:

    A = B / 2;
    */
    Tensor result(dim1, dim2, 0.0);
    int i, j;
    for (i = 0; i < this->dim1; i++)
    {
        for (j = 0; j < this->dim2; j++)
        {
            result(i, j) = this->self[i][j] / scalar;
        }
    }
    return result;
}

// Take any given matrices transpose and returns another Tensor
template <typename T>
Tensor<T> Tensor<T>::transpose()
{
    /* 
    It takes the transpose of the given Tensor.

    Example:

    A = B.transpose();
    */
    Tensor<T> Transpose(dim2, dim1, 0.0);
    for (int i = 0; i < this->dim2; i++)
    {
        for (int j = 0; j < this->dim1; j++)
        {

            Transpose(i, j) = this->self[j][i];
        }
    }
    return Transpose;
}
template <typename T>
Tensor<T> Tensor<T>::inner(Tensor<T> &in)
{
    Tensor<T> result(this->dim1, this->dim2, 0.0);
    for (int i = 0; i < this->dim1; i++)
    {
        for (int j = 0; j < this->dim1; j++)
        {
            result(i, j) = this->self[i][j] * in(i, j);
        }
    }
    return result;
}
template <typename T>
T &Tensor<T>::operator()(const int &dim1No, const int &dim2No)
{
    /*
    This operator returns the corresponding element of a given (dim1,dim2umn) tuple.

    Example:
    double a;
    Tensor<double> A(5,5,5.0);

    a = A(2,4);
    */
    return this->self[dim1No][dim2No];
}
template <typename T>
const T &Tensor<T>::operator()(const int &dim1No, const int &dim2No) const
{
    /*
    This operator returns the corresponding element of a given (dim1,dim2umn) tuple.

    Example:
    double a;
    Tensor<double> A(5,5,5.0);

    a = A(2,4);
    */
    return this->self[dim1No][dim2No];
}

template <typename T>
T &Tensor<T>::operator()(const int &dim1No)
{
    /*
    This operator returns the corresponding element of a given (dim1) number of a std::vector.

    Example:
    double a;
    Tensor<double> A(5,5.0);

    a = A(2);
    */
    if (this->dim2 == 1)
    {
        return this->self[dim1No][0];
    }
}
template <typename T>
const T &Tensor<T>::operator()(const int &dim1No) const
{
    /*
    This operator returns the corresponding element of a given (dim1) number of a std::vector.

    Example:
    double a;
    Tensor<double> A(5,5.0);

    a = A(2);
    */
    if (this->dim2 == 1)
    {
        return this->self[dim1No][0];
    }
}

template <typename T>
Tensor<int> Tensor<T>::getSize() const
{
}

template <typename T>
std::vector<std::vector<T>> Tensor<T>::getData() const
{
    /*
    This function returns self 2d std::vector.

    Example:
    std::vector<std::vector<double>> a;
    Tensor<double> A(5,5,5.0);

    a = A.getMat();
    */
    return this->self;
}

// Prints the Tensor beautifully
template <typename T>
void Tensor<T>::print() const
{
    /*
    This function prints the Tensor in a beautiful and understanble way.

    Example:
    
    Tensor<int> A(2,2,3);

    A.print() =====> [3] [3]
                     [3] [3]
    */
    std::cout << "Tensor: " << endl;
    for (int i = 0; i < this->dim1; i++)
    {
        for (int j = 0; j < this->dim2; j++)
        {
            std::cout << "[" << this->self[i][j] << "] ";
        }
        std::cout << endl;
    }
}