#include <iostream>
#include <cmath>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"

using namespace std;

xt::xarray<double> spa_decode(xt::xarray<double> r, xt::xarray<int> H){
    bool stop = false;
    int Imax = 10;
    int I = 0;
    auto sh = H.shape();
    xt::xarray<double> M = xt::zeros<double>(H.shape());
    xt::xarray<double> E = xt::zeros<double>(H.shape());
    xt::xarray<double> l = xt::zeros<double>(r.shape());
    xt::xarray<int> H_mirr = (H + xt::ones<int> (H.shape())) % 2;

    while (stop == false and I != Imax){
        if (I == 0){
            for (int j = 0; j < sh[0]; j++){
            xt::row(M, j) = r*xt::row(H, j);
            }
        }
        cout << "M:" << endl << M << endl;

        //calc E
        xt::xarray<double> M1 = xt::tanh(M/2) + H_mirr;
        for (int j = 0; j < sh[0]; j++){
            for (int i = 0; i < sh[1]; i++){
                if (H(j, i) != 0){
                    xt::xarray<double> res = xt::prod(xt::row(M1, j)) / M1(j, i);
                    E(j, i) = log((1 + res(0)) / (1 - res(0)));
                }
            }
        }
        cout << "E:" << endl << E << endl;

        l = r + xt::sum(E, {0});
        cout << "l:" << endl << l << endl;
        // NRZ
        for (int i = 0; i < l.size(); i++){
            if (l(i) > 0) {
                l(i) = 0;
            }
            else l(i) = 1;
        }
        cout << "decoded:" << endl << l << endl;

        xt::xarray<int> s = xt::linalg::dot(H, l);
        s = s % 2;
        cout << "s:" << s << endl;

        if (s == xt::zeros<int>(s.shape())) {
            stop = true;
        }
        else {
            I++;
            //calc M
            for (int j = 0; j < sh[0]; j++){
                for (int i = 0; i < sh[1]; i++){
                    if (H(j, i) != 0){
                        auto tmp = xt::sum(xt::col(E, i));
                        M(j, i) = tmp[0] - E(j, i) + r(i);
                    }
                }
            }
            M = M * H;
        }
    }
  return l;
}

int main(){
    xt::xarray<double> H
    {{1, 1, 0, 1, 0, 0},
    {0, 1, 1, 0, 1, 0},
    {1, 0, 0, 0, 1, 1},
    {0, 0, 1, 1, 0, 1}};

    xt::xarray<double> r1 // [0 0 1 0 1 1]
    {-1.3863, 1.3863, -1.3863, 1.3863, -1.3863, -1.3863};

    xt::xarray<double> r2 // [0 0 1 0 1 1]
    {-0.5, 2.5, -4.0, 5.0, -3.5, 2.5};

    xt::xarray<double> l = spa_decode(r2, H);
    cout << "Decoded Message: " << l << endl;
    return 0;
}

