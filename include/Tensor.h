#ifndef TENSOR
#define TENSOR
/* 
It is the template class header file of the Tensor class. This class has a lot of different constructor schemas. The template nature of this
class allows us to work with different data types. For example, in training process, we can work with the double and integer data types, while, in
sampling process, we can work with the chars or char arrays. Also, it can be used in ordered maps. In this way, we can construct the whole
system as a closed box. 
*/

// TODO: Template Specialization for the char and int data types. It may be useful for future works.

#include <vector>
#include <string>

template <typename T>
class Tensor
{

public:
	Tensor<T>();
	Tensor<T>(const int, const int, const int, const std::string);
	Tensor<T>(const int, const int, const std::string);
	Tensor<T>(const int, const std::string);
	Tensor<T>(const Tensor<T> &);
	~Tensor<T>();

	Tensor<T> operator=(const Tensor<T> &);
	Tensor<T> operator+(const Tensor<T> &);
	Tensor<T> operator+=(const Tensor<T> &);
	Tensor<T> operator-(const Tensor<T> &);
	Tensor<T> operator-=(const Tensor<T> &);
	Tensor<T> operator*(const Tensor<T> &);
	Tensor<T> operator*=(const Tensor<T> &);

	Tensor<T> operator+(T);
	Tensor<T> operator-(T);
	Tensor<T> operator*(T);
	Tensor<T> operator/(T);

	Tensor<T> transpose();
	Tensor<T> inner(Tensor<T> &);
	Tensor<T> concat(const Tensor<T> &, const int);

	T &operator()(const int &, const int &);
	const T &operator()(const int &, const int &) const;
	T &operator()(const int &);
	const T &operator()(const int &) const;

	Tensor<int> getSize() const;
	std::vector<std::vector<T>> getData() const;

	void print() const;

private:
	int dim1;
	int dim2;
	int dim3;
	std::vector<std::vector<std::vector<T>>> self;
};

#include "../source/Tensor.cpp"

#endif