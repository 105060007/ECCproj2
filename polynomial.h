#include <iostream>
#include <vector>

const int MAX = 32;     // degree: 0~31

class Polynomial{
public:
	Polynomial();
	~Polynomial();
	void setPolynomial(int coefficient[]);
	void printPolynomial();
    int getDegree();
	std::vector<int> getCoeff();

	Polynomial operator+(const Polynomial addi);
	Polynomial operator*(const Polynomial mult);

private:
	int coeff[MAX];     // start w/ low degree
	int degree;
};