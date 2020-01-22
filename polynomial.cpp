#include "polynomial.h"
#include "alpha.h"
using namespace std;

Polynomial::Polynomial(){
	for (int i=0; i < MAX; i++)
		coeff[i] = 0;
    degree = 0;
}

Polynomial::~Polynomial(){
}

void Polynomial::setPolynomial(int coefficient[]){
	for(int i=0; i < MAX; i++)
        this->coeff[i] = coefficient[i];
    
    // search nonzero coeff from high degree
	for (int i=MAX-1; i >=0; i--){
		if(coeff[i] != 0){
            this->degree = i;
            break;
        }
    }
}

void Polynomial::printPolynomial(){
    if(degree == 0) cout << coeff[0] << endl;
    else{
        if(coeff[0] != 0) cout << coeff[0] << " + ";
        for (int i=1; i < degree; i++)
            if(coeff[i]!=0) cout << coeff[i] << "x^" << i << " + ";
        cout << coeff[degree] << "x^" << degree << endl;
    }
}

int Polynomial::getDegree(){
    return degree;
}

vector<int> Polynomial::getCoeff(){
    vector<int> coefficient;
    for(int i=0; i <= degree; i++)
        coefficient.push_back(coeff[i]);
    return coefficient;
}

Polynomial Polynomial::operator+(const Polynomial addi){
    Polynomial poly;
    int tmp[32] = {0};

	for (int i=0; i <= max(this->degree, addi.degree); i++){
		tmp[i] = (this->coeff[i] ^ addi.coeff[i]);
	}
    poly.setPolynomial(tmp);

	return poly;
}

Polynomial Polynomial::operator*(const Polynomial mult){
	Polynomial poly;
    int m = this->degree;
    int n = mult.degree;
    if(m+n > 31) cout << "error: degree sum > 31." << endl;

    int prod[32] = {0}; 
    int tmp;

    for (int i=0; i<=m; i++){
        for (int j=0; j<=n; j++){
            if(this->coeff[i]==0 || mult.coeff[j]==0){
                continue;
            }
            else{
                tmp = invAlpha[this->coeff[i]] + invAlpha[mult.coeff[j]];   // 算出是alpha幾次方
                prod[i+j] ^= alpha[(tmp+31)%31];            // XOR
            }
        }
    }
    poly.setPolynomial(prod);

	return poly;
}