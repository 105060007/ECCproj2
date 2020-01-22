#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <tuple>
#include <math.h>

#include "polynomial.h"
#include "alpha.h"

using namespace std;

Polynomial LongDivision(Polynomial Nume, Polynomial Deno);
tuple<Polynomial, Polynomial> Euclidean(Polynomial ax, Polynomial bx, int myu, int nyu);
int Substitute(Polynomial poly, int i);

int main(int argc, char* argv[]){
    // read file
    ifstream file(argv[1]); // pass file name as argment
	string linebuffer;
    int codeword[31] = {0};       // * = 0
    ofstream out_data("output.txt");

	while (file && getline(file, linebuffer)){
        int decodedCodeword[31];
        bool decode = true;
        vector<int> erasureSet;
		if (linebuffer.length() == 0) continue;        

		// read the codeword
        istringstream iss(linebuffer);
        string buf;
        for(int i=0; i<31; i++){
            iss >> buf;
            if(buf == "*"){
                erasureSet.push_back(i);
                codeword[i] = 0;
            }
            else codeword[i] = stoi(buf);
        }
        
        // check the # of erasure
        cout << endl << "Size of erasureSet: " << erasureSet.size() << endl;
        if(erasureSet.size() > 16){
            cout << "fail: Out of Error-Correcting range. (# of erasure = " << erasureSet.size() << ")" << endl;
            // decode = false;
            out_data << "fail: Out of Error-Correcting range. (# of erasure = " << erasureSet.size() << ")" << endl;
            continue;
        }

        // decode!
        cout << "Start decodeing..." << endl;
        // calculate the sydrome
        int syndrome[32] = {0};
        for(int j=1; j<=16; j++){
            for(int i=0; i<31; i++){
                if(codeword[i]==0) continue;
                else{
                    int exp = (invAlpha[codeword[i]]+i*j)%31;
                    syndrome[j-1] = syndrome[j-1] ^ alpha[exp];
                }
            }
        }
        Polynomial Syndrome;
        Syndrome.setPolynomial(syndrome);
        cout << "Syndrome: ";
        Syndrome.printPolynomial();

        // calculate the erasure locater poly
        Polynomial erasureLocator;
        Polynomial tmpPoly;

        if(erasureSet.empty()){
            int a[32] = {1, 0};
            erasureLocator.setPolynomial(a);
        }
        else{
            int a[32] = {1, alpha[erasureSet[0]]};
            erasureLocator.setPolynomial(a);
            for(int i=1; i<erasureSet.size(); i++){
                a[1] = alpha[erasureSet[i]];
                tmpPoly.setPolynomial(a);
                erasureLocator = erasureLocator * tmpPoly;
            }
        }
        cout << "erasureLocator: ";
        erasureLocator.printPolynomial();

        // calculate S0(x) and x^r, then start Euclid's algo.
        Polynomial S0 = erasureLocator * Syndrome;
        int tmp[32] = {0}; tmp[16] = 1;
        Polynomial xr; xr.setPolynomial(tmp);

        tuple<Polynomial, Polynomial> err;
        int myu = (16-erasureSet.size())/2;
        int nyu = ceil((16.0+erasureSet.size())/2) -1;
        if(S0.getDegree() >= 16){
            vector<int> modXr = S0.getCoeff();
            int modXR[32] = {0};
            for(int i=0; i<16; i++) modXR[i] = modXr[i];
            S0.setPolynomial(modXR);
        }
        
        cout << "S0: " ;
        S0.printPolynomial();
        cout << "xr: ";
        xr.printPolynomial();
        err = Euclidean(xr, S0, myu, nyu);

        Polynomial errorLocator, errorEval; 
        errorLocator = get<0>(err);     // sigma1(x) 
        errorEval = get<1>(err);        // omega(x)
        cout << "errorLocator: ";
        errorLocator.printPolynomial();
        cout << "errorEval: ";
        errorEval.printPolynomial();

        // Time domain Approach
        Polynomial sigma = erasureLocator * errorLocator;
        Polynomial sigmaDerivative;
        int E[31] = {0};                // R+E = R'
        // formal derivative?
        vector<int> d = sigma.getCoeff();
        tmp[16] = 0;
        for(int i=1; i<d.size(); i++){
            tmp[i-1] = (i%2==0) ? 0 : d[i];
        }
        sigmaDerivative.setPolynomial(tmp);
        cout << "sigma: ";
        sigma.printPolynomial();
        // sigmaDerivative.printPolynomial();

        bool check = false;
        if(d.front()==1 && d.size()==1) check = true;
        if(check){
            // no error and erasure
            for(int i=0; i<31; i++){
                cout << codeword[i] << " ";
                out_data << codeword[i] << " ";
            }
            cout << endl;
            out_data << endl;
            continue;
        }

        int count = 0;
        if(d.front()==0){
            decode = false;
            cout << "fail: sigma0 = 0" << endl;
            out_data << "fail: sigma0 = 0" << endl;
            continue;
        }
        else if(errorEval.getDegree()>=erasureSet.size()+sigma.getDegree()){
            decode = false;
            cout << "fail: deg(w(x)) >= e0 + deg(sigma(x))" << endl;
            out_data << "fail: der(w(x)) >= e0 + deg(sigma(x))" << endl;
            continue;
        }
        else{
            count = 0;
            for(int i=0; i<31; i++){
                int k = Substitute(sigmaDerivative, -i);
                if(Substitute(sigma, -i)==0 && k!=0){
                    count ++;
                    int w = Substitute(errorEval, -i);
                    if(w==0) E[i] = 0;
                    // E[i] = -w/k;
                    else E[i] = alpha[(invAlpha[w]-invAlpha[k]+31)%31];
                }
                else E[i] = 0;
            }

            if(count == sigma.getDegree()){
                decode = true;
                for(int i=0; i<31; i++){
                    // R = R' ^ E
                    decodedCodeword[i] = codeword[i] ^ E[i];
                    cout << decodedCodeword[i] << " ";
                    out_data << decodedCodeword[i] << " ";
                }
                cout << endl;
                out_data << endl;
            }
            else{
                decode = false;
                cout << "fail: count != deg(sigma(x))" << endl;
                out_data << "fail: count != deg(sigma(x))" << endl;
            }
        }
    }
    out_data.close();
	return 0;
}


tuple<Polynomial, Polynomial> Euclidean(Polynomial ax, Polynomial bx, int myu, int nyu){
    vector<Polynomial> u, v, r, q;
    int i;

    Polynomial tmp1, tmp2;
    int tmp[32] = {0};
    tmp1.setPolynomial(tmp);    // 0
    tmp[0] = 1;
    tmp2.setPolynomial(tmp);    // 1

    u.push_back(tmp2); u.push_back(tmp1);
    v.push_back(tmp1); v.push_back(tmp2);
    r.push_back(ax); r.push_back(bx);
    q.assign(2, tmp1);

    if(v[1].getDegree()<=myu && r[1].getDegree()<=nyu){
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
            tmp1 = q[i] * v[i-1]; v.push_back(tmp1 + v[i-2]);
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
    
	r.resize(dr+1);

    if(dN >= dD){
        while (dN >= dD){
            int d[32] = {0};
            // 將除式右移到和被除式degree一樣
			for(int i=0; i<=dD; i++) d[i+dN-dD] = D[i];
			dd = dN;
            // 兩個最高項係數相除得到商
            if(N[dN] == 0){
                q[dN-dD] = 0;
            }
            else {
                int tmp =  invAlpha[N[dN]] + 31 - invAlpha[d[dN]];
                q[dN-dD] = alpha[tmp%31];
            }
            // d[i] = d[i] * q[dN-dD]; 
            for(int i=0; i<32; i++){
                if(d[i]==0 || q[dN-dD]==0) {d[i] = 0;}
                else{
                    int tmp = invAlpha[d[i]] + invAlpha[q[dN-dD]];
                    d[i] = alpha[tmp%31];
                }
            }
            // 被除式 減掉 除式*商，degree -1
            for(int i=0 ; i<=N.size(); i++ ) {N[i] = N[i] ^ d[i];}
            dN--;
        }
        
        Polynomial quo;
        quo.setPolynomial(q);

        return quo;
    }
    else{
        cout << "fail: degree of Denominator > degree of Numerator." << endl;
        return Deno;
    }
}

int Substitute(Polynomial poly, int i){
    vector<int> coeff = poly.getCoeff();
    int result = 0;
    i = (i+31)%31;
    for(int j=0; j<=poly.getDegree(); j++){
        // coeff[i] * (alpha[-i])^j
        int tmp;
        if(coeff[j]!=0){
            tmp = (invAlpha[coeff[j]] + i*j)%31;
            result = result ^ alpha[tmp];
        }
    }

    return result;
}