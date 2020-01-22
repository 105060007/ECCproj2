#include "polynomial.h"
#include "alpha.h"
#include <tuple>

using namespace std;
Polynomial LongDivision(Polynomial Nume, Polynomial Deno);
tuple<Polynomial, Polynomial> Euclidean(Polynomial ax, Polynomial bx, int myu, int nyu);

int main(){
	// Polynomial ax, bx;
	// Polynomial ErrorLocator, ErrorEval;

    // int a[32] = {0}; a[8] = 1;
    // int b[32] = {0}; b[0]=1;b[1]=1;b[2]=1;b[4]=1;b[6]=1;
    // ax.setPolynomial(a);
    // bx.setPolynomial(b);

	// tuple<Polynomial, Polynomial> err;
	// err = Euclidean(ax, bx, 5, 2);

	// ErrorLocator = get<0>(err);
	// ErrorEval = get<1>(err);
	// cout << "Outcome:" << endl;
	// ErrorLocator.printPolynomial();
	// ErrorEval.printPolynomial();

    Polynomial quo;
    Polynomial S0;
    int tmpS0[32] = {21, 24, 4, 18, 16, 9, 16, 27, 28, 30, 30, 25, 25, 17, 13, 31};
    S0.setPolynomial(tmpS0);
    int tmp[32] = {0}; tmp[16] = 1;
    Polynomial xr; xr.setPolynomial(tmp);
    
    quo = LongDivision(xr, S0);
    cout << "Outcome: ";
    quo.printPolynomial();
}

tuple<Polynomial, Polynomial> Euclidean(Polynomial ax, Polynomial bx, int myu, int nyu){
    vector<Polynomial> u, v, r, q;

    Polynomial tmp1, tmp2;
    int tmp[32] = {0};
    tmp1.setPolynomial(tmp);    // 0
    tmp[0] = 1;
    tmp2.setPolynomial(tmp);    // 1

    u.push_back(tmp2); u.push_back(tmp1);
    v.push_back(tmp1); v.push_back(tmp2);
    r.push_back(ax); r.push_back(bx);
    q.assign(2, tmp1);

    int i;

    // 到底要檢查啥
    if(v[1].getDegree()<=myu && r[1].getDegree()<=nyu){
        cout << "if~~~???" << endl;
        tuple<Polynomial, Polynomial> myTuple;
        get<0>(myTuple) = v[1];
        get<1>(myTuple) = r[1];
        
        return myTuple;
    }
    else{
        i = 1;
        do{
            i++;
            // calculate q[i]
            q.push_back(LongDivision(r[i-2], r[i-1]));

            // update
            tmp1 = q[i] * r[i-1]; r.push_back(tmp1 + r[i-2]);
            tmp1 = q[i] * u[i-1]; u.push_back(tmp1 + u[i-2]);
            tmp1 = q[i] * v[i-1];  v.push_back(tmp1 + v[i-2]);
        } while ((v[i].getDegree() > myu) || (r[i].getDegree() > nyu));

        tuple<Polynomial, Polynomial> myTuple;
        get<0>(myTuple) = v[i];
        get<1>(myTuple) = r[i];
        
        return myTuple;
    }
}

Polynomial LongDivision(Polynomial Nume, Polynomial Deno){

    // N = D*q + r
    int dN, dD, dd, dq, dr;
    vector<int> N = Nume.getCoeff();
    vector<int> D = Deno.getCoeff();

    for(int i=0; i<=16; i++) cout << N[i] << " ";
    cout << endl;
    for(int i=0; i<=15; i++) cout << D[i] << " ";
    // vector<int> q, r, d;
    vector<int> r;
    int q[32] = {0};
    dN = N.size() -1;
    dD = D.size() -1;
    dq = dN-dD;
    dr = dN-1;

    if( dN < 1 || dN < 1 ) {
		cerr << "Error: degree of D and N must be positive.\n";
		return Nume;
	}
    
	// d.resize(dN+1);
    // q.resize(dq+1);
	r.resize(dr+1);

    if(dN >= dD){
        while (dN >= dD){
            // d.assign(d.size(), 0);
            int d[32] = {0};
 
			for(int i = 0 ; i <= dD ; i++ )
				d[i+dN-dD] = D[i];
			dd = dN;

            if(N[dN] == 0){
                q[dN-dD] = 0;
            }
            else {
                int tmp =  invAlpha[N[dN]] + 31 - invAlpha[d[dd]];
                q[dN-dD] = alpha[tmp%31];
                cout << invAlpha[N[dN]] << " " << invAlpha[d[dd]] << endl;
            }
            cout << "N[dN] = " << N[dN] << ", d[dd] = " << d[dd] << endl; 
            cout << "q = " << q[dN-dD] << endl;

            for(int i=0 ; i<32; i++){
                // d[i] = d[i] * q[dN-dD];         // alpha?
                if(d[i]==0 || q[dN-dD]==0) {d[i] = 0;}
                else{
                    int tmp = invAlpha[d[i]] + invAlpha[q[dN-dD]];
                    d[i] = alpha[tmp%31];
                }
            }
                           
            for(int i=0 ; i<=N.size(); i++ )
                N[i] = N[i] ^ d[i];             // alpha?
            dN--;
        }
    }
    else{
        cout << "fail: degree of Denominator > degree of Numerator." << endl;
        return Nume;
    }

    // int qCoeff[32] = {0};
    // Polynomial quo;
    // for(int i=0; i<q.size(); i++) qCoeff[i] = q[i];
    // quo.setPolynomial(qCoeff);
    Polynomial quo;
    quo.setPolynomial(q);

    return quo;
}