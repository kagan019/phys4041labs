#include <bits/stdc++.h>
using namespace std;

int ipow(int b, int e) {
    int s = 1;
    for (int i = 0; i < e; i++)
        s *= b;
    return s;
}

//a)

float form1(int n) {
    return ipow(-1,n)*((float)n)/(n+1);
}

float form2a(int n) {
    return (float)(2*n-1) / -(2*n);
}
float form2b(int n) {
    return (float)(2*n) / (2*n+1);
}

float form3(int n) {
    return 1.0f/(2*n*(2*n+1));
}

int main() {
    cout << fixed << setprecision(15);
    vector<int> Ns = {100,1000,10000,2500000};
    for (int &N : Ns) { 
        //for each trial,
        cout << "N = " << N << "\n";
        float s1 = 0,s2a = 0, s2b =0, s3 = 0;
        for (int n = 1; n < 2*N;n++) {
            s1 += form1(n);
        }
        for (int n = 1; n < N;n++) {
            s2a += form2a(n);
            s2b += form2b(n);
            s3 += form3(n);
        }
        float s2 = s2a+s2b;
        cout << "\ts1 = " << s1 << endl; 
        cout << "\ts2 = " << s2 << endl; 
        cout << "\ts3 = " << s3 << endl; 
    }
    return 0;
}


/* comments

s1 across trials is accurate to 2 decimal places. but it is grossly inaccurate to s2 and s3,
  suggesting that this method is not numerically stable.

s2 for N=100 is accurate to 2 decimal places with s2 for N=1000 and N=10000, 
  and those two are accurate to eachother for all 15 decimal places (which is 
  well beyond the precision of singles). This series appears to be converging 
  quickly
  
s3 and s2 for N = 100 are accurate to 4 decimal places. s3 for N = 100 is accurate 
  to s3 for N = 1000 also to 2 decimal places. s3 for N=1000 and 10000 compare similar 
  only to 3 decimal places. Alarmingly, s3 for N=10000 has overshot the value s2 appears
  to be converging to by 0.002

*/

//b)
// we need to know how many terms k until the series from n=k to infinity is less than 10^-7
// we approximate this series with an integral over continuous variable x, for an overestimate
// and thus an upper bound on error. We can do this because the formula is positive and monotonically
// decreasing

/*
series from k to infinity of 1/(2n(2n+1)) <= integral from k to infinity of dx/(2x(2x+1)) <= 10^-7
PFD
1/(2x(2x+1)) = A/2x + B/(2x+1)
1 = A(2x+1) + B(2x)
clearly,
A = 1 
B = -1

1/2 * int_k^inf dx/x - 1/2 * int_k^inf du/u 
                            du = 2dx  u = 2x+1
=1/2 (ln|2x| - ln|2x+1|)|_k^inf
first:
1/2 (ln(2x/(2x+1)) as x->inf = 0
    2x/(2x+1) clearly approaches 1 as x->inf (by L'Hopital's rule)
so we just need to find
ln((2k+1)/(2k)) = 2 * 10^-7
2k+1 = e^(...) * 2k
1 = 2k(e^(..)-1)
k = 1/(2(e^(2*10^-7) - 1)) ~= 2500000
*/